"""
Test dataset comparing grism analysis codes (aXe, Pirzkal, Brammer)

http://archive.stsci.edu/hst/search.php?action=Search&sci_instrume=WFC3&sci_data_set_name=IBHM52*
"""

def make_catalog():
    """
    Make F160W-selected catalog from the CANDELS reference image
    """
    import pyraf
    from pyraf import iraf
    
    ######
    ## Make smaller mosaic image
    # swarp /Users/brammer/3DHST/Ancillary/Mosaics/cosmos_3dhst.v4.0.F160W_orig_sci.fits -c subimage.swarp 
        
    ######
    ## Run SExtractor
    ROOT = 'cosmos_F160W'
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['CHECKIMAGE_NAME'] = '%s_seg.fits' %(ROOT)
    se.options['WEIGHT_TYPE']     = 'NONE'
    se.options['GAIN'] = '0'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_THRESH']    = '1.8' 
    se.options['ANALYSIS_THRESH']  = '1.8' 
    se.options['DETECT_MINAREA'] = '15' 
    se.options['MASK_TYPE'] = 'NONE'
    se.options['DEBLEND_NTHRESH']    = '32' 
    se.options['DEBLEND_MINCONT']  = '0.001' 
    se.options['SEEING_FWHM'] = '0.12'
    
    se.options['BACK_TYPE'] = 'AUTO'
    se.options['BACKPHOTO_TYPE'] = 'LOCAL'
    
    se.options['GAIN'] = '3000.'

    se.options['MAG_ZEROPOINT'] = '25.96'
    se.options['CATALOG_TYPE']    = 'ASCII_HEAD'
    se.options['CATALOG_NAME']    = '%s.cat' %(ROOT)
    status = se.sextractImage('cosmos_F160W_mosaic_sci.fits[0]')
    threedhst.sex.sexcatRegions('%s.cat' %(ROOT), '%s.reg' %(ROOT), format=2)
    
def fit_separate():
    """
    Make separate spectra for each FLT
    """
    sub = 'abcd'
    sf = threedhst.shifts.ShiftFile('cosmos-24-F140W_shifts.txt')
    
    #### Preparation steps
    for i in range(4):
        for filt in ['F140W','G141']:
            asn = threedhst.utils.ASNFile('cosmos-24-%s_asn.fits' %(filt))
            asn.product = asn.product.replace('24','24%s' %(sub[i])).lower()
            asn.exposures = [asn.exposures[i]]
            print asn.product
            asn_file = asn.product.lower() + '_asn.fits'
            asn.write(asn_file)
            ### Shift file
            threedhst.shifts.make_blank_shiftfile(asn_file, xshift=sf.xshift[0], yshift=sf.yshift[0] ,rot=sf.rotate[0])
            ### Distortion corrected
            threedhst.prep_flt_files.startMultidrizzle(asn_file, use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=1, driz_cr=False, updatewcs=False, clean=True, median=False, build_drz=True)
    
    #### "interlaced" images
    REF_ROOT='cosmos_F160W'
    NGROW=150
    for i in range(1,4):
        pointing = 'cosmos-24%s' %(sub[i])
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT=pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=False, pad=60, REF_ROOT=REF_ROOT, CATALOG='../Catalog/cosmos_F160W.cat',  NGROW=NGROW, verbose=True, auto_offsets=True, growx=1, growy=1)
        #########
        ## Generate the 1x1 interlaced image. 
        NGROW=150
        unicorn.reduce.interlace_combine('%s-G141' %(pointing), pad=60, NGROW=NGROW, growx=1, growy=1, view=False)
        unicorn.reduce.interlace_combine('%s-F140W' %(pointing), pad=60, NGROW=NGROW, growx=1, growy=1, view=False)
        #########
        ## Be sure to use the most recent configuration file
        unicorn.reduce.set_grism_config(force=True)
        #########
        ## Generate the grism model
        model = unicorn.reduce.process_GrismModel(root=pointing, grow_factor=1, growx=1, growy=1, MAG_LIMIT=24, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, model_slope=0, direct='F140W', grism='G141')
        #########
        ## Fit an object
        id = 2324    
        model.twod_spectrum(id, miny=-17, USE_REFERENCE_THUMB=True)
        gris = unicorn.interlace_fit.GrismSpectrumFit('%s_%05d' %(model.root, id), skip_photometric=False)
        gris.fit_in_steps(zrfirst=[1.2,1.4], oned=False)
    
    import time
    times = []
    models = []
    for i in range(4):
        pointing = 'cosmos-24%s' %(sub[i])
        t0 = time.time()
        model = unicorn.reduce.process_GrismModel(root=pointing, grow_factor=1, growx=1, growy=1, MAG_LIMIT=26, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, model_slope=0, direct='F140W', grism='G141')
        t1 = time.time()
        times.append(t1-t0)
        models.append(model)
        #
        model.refine_mask_background(update=True, threshold=0.004, grow_mask=12)
    
    #
    id = 1260
    for model in models:
        model.twod_spectrum(id, miny=-34, USE_REFERENCE_THUMB=True)
        gris = unicorn.interlace_fit.GrismSpectrumFit('%s_%05d' %(model.root, id), skip_photometric=False)
        zr = np.clip([gris.z_peak-0.1, gris.z_peak+0.1], 0,4)
        gris.fit_in_steps(zrfirst=zr, oned=True)
       
    #
    full_model = unicorn.reduce.process_GrismModel(root='cosmos-24', grow_factor=2, growx=2, growy=2, MAG_LIMIT=24, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, model_slope=0, direct='F140W', grism='G141')
    full_model.twod_spectrum(id, miny=-34, USE_REFERENCE_THUMB=True)
    
    full = unicorn.interlace_fit.GrismSpectrumFit('%s_%05d' %('cosmos-24', id), skip_photometric=False)
    full.fit_in_steps(zrfirst=zr, oned=True)
    
    # %timeit gris.twod.compute_model()
    # %timeit full.twod.compute_model()
    
def prep():
    
    import unicorn
    
    ########
    ## Initial preparation: align to reference image and subtract backgrounds
    ALIGN = '/Users/brammer/3DHST/Ancillary/Mosaics/cosmos_3dhst.v4.0.F160W_orig_sci.fits'
    
    unicorn.candels.make_asn_files(uniquename=False)
    
    asn_files = glob.glob('*F140W_asn.fits')
    for file in asn_files:
        if os.path.exists(file.replace('asn','drz')):
            continue
        ## First preparation steps
        threedhst.prep_flt_files.process_3dhst_pair(file, file.replace('F140W', 'G141'), adjust_targname=False, ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate, shift')
        
    #########
    ## Blot reference image and catalog to interlaced frame
    REF_ROOT='cosmos_F160W'
    unicorn.reduce.prepare_blot_reference(REF_ROOT=REF_ROOT, filter='F160W', REFERENCE='../Catalog/cosmos_F160W_mosaic_sci.fits', SEGM = '../Catalog/cosmos_F160W_seg.fits', sci_extension=0)

    pointing = 'cosmos-24d'
    unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT=pointing+'-F140W', NGROW=NGROW, verbose=True)
    unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=False, pad=60, REF_ROOT=REF_ROOT, CATALOG='../Catalog/cosmos_F160W.cat',  NGROW=NGROW, verbose=True, auto_offsets=True)
    
    #########
    ## Generate the 2x2 interlaced image.  The NGROW parameter adds some space
    ## on the edge of the FLT frame to get objects not in the direct image
    ## but nearby on the reference image
    NGROW=150
    unicorn.reduce.interlace_combine('%s-G141' %(pointing), pad=60, NGROW=NGROW, growx=2, growy=2, view=False)
    unicorn.reduce.interlace_combine('%s-F140W' %(pointing), pad=60, NGROW=NGROW, growx=2, growy=2, view=False)
    
    #########
    ## Be sure to use the most recent configuration file
    unicorn.reduce.set_grism_config(force=True)
    
    #########
    ## Generate the grism model
    model = unicorn.reduce.process_GrismModel(root='cosmos-24', grow_factor=2, growx=2, growy=2, MAG_LIMIT=26, REFINE_MAG_LIMIT=21.5, make_zeroth_model=False, use_segm=False, model_slope=0, direct='F140W', grism='G141')
    
    im = pyfits.open('cosmos-24_inter_model.fits')
    img = pyfits.open('cosmos-24-G141_inter.fits')
    asn = threedhst.utils.ASNFile('cosmos-24-G141_asn.fits')
    
    flt_0 = img[1].data[NGROW*2+30:-(NGROW*2+30):2, NGROW*2+30:-(NGROW*2+30):2]*4
    model_0 = im[0].data[NGROW*2+30:-(NGROW*2+30):2, NGROW*2+30:-(NGROW*2+30):2]*4
    
    pyfits.writeto('ibhm52c4q_flt_CONT_2x2.fits', data=model_0, clobber=True)
    
    ## check background subtraction
    model.refine_mask_background(update=True, threshold=0.001, grow_mask=30)
    
    t0 = time.time()
    model = unicorn.reduce.process_GrismModel(root='cosmos-24a', grow_factor=1, growx=1, growy=1, MAG_LIMIT=24, REFINE_MAG_LIMIT=21.5, make_zeroth_model=False, use_segm=False, model_slope=0, direct='F140W', grism='G141')
    t1 = time.time()
    print (t1-t0)/316  # (316 total H < 24, 81 refined H<21.5)  127.2 s total, 0.4 s / object
    
    #########
    ## Extract a spectrum and fit it
    # id = 1260
    # id = 1076
    # id = 1348
    id = 2324    
    model.twod_spectrum(id, miny=-34, USE_REFERENCE_THUMB=True)
    gris = unicorn.interlace_fit.GrismSpectrumFit('%s_%05d' %(model.root, id), skip_photometric=False)
    gris.fit_in_steps(zrfirst=[1.2,1.4], oned=True)
    
def strip_flt():
    """
    Pull model FLTs out of padded model.fits images
    """
    NGROW, pad = 150, 60
    off = NGROW+pad/2
    for sub in 'abcd':
        root='cosmos-24%s' %(sub)
        im = pyfits.open('%s_inter_model.fits' %(root))
        asn = threedhst.utils.ASNFile('%s-G141_asn.fits' %(root))
        flt_im = im[0].data[off:-off, off:-off]
        flt = pyfits.open('%s_flt.fits' %(asn.exposures[0]))
        pyfits.writeto('%s_flt_CONT.fits' %(asn.exposures[0]), data=flt_im, header=flt[1].header, clobber=True)
    
def compare_spectra():
    import mywfc3.stgrism as st
    import unicorn
    
    ### Fancy colors
    import seaborn as sns
    import matplotlib.pyplot as plt
    cmap = sns.cubehelix_palette(as_cmap=True, light=0.95, start=0.5, hue=0.4, rot=-0.7, reverse=True)
    cmap.name = 'sns_rot'
    plt.register_cmap(cmap=cmap)
    sns.set_style("ticks", {"ytick.major.size":3, "xtick.major.size":3})
    plt.set_cmap('sns_rot')
    #plt.gray()
    
    fig = st.compare_methods(x0=787, y0=712, v=np.array([-1.5,4])*0.6, NX=180, NY=40, direct_off=100, final=True, mask_lim = 0.02)
    #fig.tight_layout()
    unicorn.plotting.savefig(fig, '/tmp/compare_model_star.pdf', dpi=300)

    fig = st.compare_methods(x0=485, y0=332, v=np.array([-1.5,4])*0.2, NX=180, NY=40, direct_off=100, final=True, mask_lim = 0.1)
    unicorn.plotting.savefig(fig, '/tmp/compare_model_galaxy.pdf', dpi=300)

    fig = st.compare_methods(x0=286, y0=408, v=np.array([-1.5,4])*0.08, NX=180, NY=40, direct_off=100, final=True, mask_lim = 0.1)
    unicorn.plotting.savefig(fig, '/tmp/compare_model_galaxy2.pdf', dpi=300)

    fig = st.compare_methods(x0=922, y0=564, v=np.array([-1.5,4])*0.2, NX=180, NY=40, direct_off=100, final=True, mask_lim = 0.15)
    unicorn.plotting.savefig(fig, '/tmp/compare_model_galaxy3.pdf', dpi=300)
    

def compare_methods(x0=787, y0=712, NX=180, NY=40, direct_off=100, final=False, v=[-1.5,4], mask_lim = 0.02):
    """
    Make a figure comparing the grism models
    """
    import pyfits
    import matplotlib.pyplot as plt
    import numpy as np

    import unicorn
    
    flt = pyfits.open('PREP_FLT/ibhm52c4q_flt.fits')
    flt_d = pyfits.open('PREP_FLT/ibhm52c2q_flt.fits')
    
    bad = (flt_d['DQ'].data & (4+4096)) > 0
    flt_d['SCI'].data[bad] = 0
    bad = (flt['DQ'].data & (4+4096)) > 0
    flt['SCI'].data[bad] = 0
         
    nor_gauss140 = pyfits.open('Nor/1501109_0023102_149_G141.GaussF140W/FLT/ibhm52c4q_flt_2.CONT.fits')
    nor_gaussAll = pyfits.open('Nor/1501109_0023102_149_G141.GaussAll/FLT/ibhm52c4q_flt_2.CONT.fits')
    nor_seg140 = pyfits.open('Nor/1501109_0023102_149_G141.SegAll/FLT/ibhm52c4q_flt_2.CONT.fits')
    nor_segAll = pyfits.open('Nor/1501109_0023102_149_G141.SegF140W/FLT/ibhm52c4q_flt_2.CONT.fits')
    
    gbb = pyfits.open('PREP_FLT/ibhm52c4q_flt_CONT_2x2.fits')
    
    sl_x_gris = slice(x0-NX/2, x0+NX/2)
    sl_x_dir = slice(x0-NX/2-direct_off, x0+NX/2-direct_off)
    sl_y = slice(y0-NY/2, y0+NY/2)
    
    direct_cutout = flt_d['SCI'].data[sl_y, sl_x_dir]
    gris_cutout = flt['SCI'].data[sl_y, sl_x_gris]
    gris_err_cutout = flt['ERR'].data[sl_y, sl_x_gris]

    models = []
    for model_im in [nor_gauss140[1].data, nor_gaussAll[1].data, nor_seg140[1].data, nor_segAll[1].data, gbb[0].data]:
        models.append(model_im[sl_y, sl_x_gris])
    
    aspect = (7*NY)/(2.*NX)
    
    fig = unicorn.plotting.plot_init(aspect=aspect, square=True, wspace=0.1, hspace=0.01, left=0.03, bottom=0.08, use_tex=final, NO_GUI=final, xs=7, fontsize=11)
    
    #v = [-1.5,4]
    lx, ly, ha, va = 0.02, 0.04, 'left', 'bottom'
    
    ax = fig.add_subplot(721)
    ax.imshow(direct_cutout, interpolation='Nearest', aspect='auto', vmin=v[0]*4, vmax=v[1]*4)
    ax.set_xticklabels([]); ax.set_yticklabels([])
    ax.text(lx, ly, 'F140W', ha=ha, va=va, transform=ax.transAxes, color='white')
    sh = direct_cutout.shape
    ax.set_xticks([0,sh[1]]); ax.set_yticks([0,sh[0]])
    
    ax.text(1-lx, 1-ly, '(%d, %d)' %(x0, y0), ha='right', va='top', transform=ax.transAxes, color='white')
    
    ax = fig.add_subplot(723)
    ax.imshow(gris_cutout, interpolation='Nearest', aspect='auto', vmin=v[0], vmax=v[1])

    mmask = (models[4] >  models[0].max()*mask_lim)
    ax.contour(mmask*1., levels=[0.9], colors='r', alpha=0.5)
    ax.set_xticklabels([]); ax.set_yticklabels([])
    ax.text(lx, ly, 'G141', ha=ha, va=va, transform=ax.transAxes, color='white')
    ax.set_xticks([0,sh[1]]); ax.set_yticks([0,sh[0]])
    
    labels = ['aXe, gauss140', 'aXe, gaussAll', 'aXe, seg140', 'aXe, segAll', 'threedhst, seg']
    
    for i in range(5):
        ax = fig.add_subplot(7,2,5+i*2)
        ax.imshow(gris_cutout-models[i], interpolation='Nearest', aspect='auto', vmin=v[0], vmax=v[1])
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(lx, ly, labels[i], ha=ha, va=va, transform=ax.transAxes, color='white')
        ax.set_xticks([0,sh[1]]); ax.set_yticks([0,sh[0]])
    
    ### Show residual histograms    
    #ax = fig.add_subplot(321)
    ax = fig.add_axes((0.52, 0.08, 0.45, 0.4))
    
    #ax = fig.subplot2grid((6, 2), (0, 1), rowspan=2)
    
    ok = (gris_cutout != 0) & (gris_err_cutout > 0) & mmask
    for i in range(5):
        ax.hist(((gris_cutout-models[i])/gris_err_cutout)[ok], bins=50, range=(-10, 10), alpha=0.8, linewidth=2, histtype='step', normed=True, label=labels[i])
    
    xx = np.arange(-10,10,0.1)
    yy = 1/np.sqrt(2*np.pi)*np.exp(-xx**2/2.)
    ax.fill_between(xx, yy, yy*0, color='black', alpha=0.2, zorder=-10)
    ax.legend(loc='upper right', prop=dict(size=8))
    ax.set_yticklabels([])
    ax.set_xlim(-8,8)
    ax.set_ylim(0,yy.max())
    ax.set_xlabel(r'Residual ($\sigma$)')
        
    #ax = fig.add_subplot(326)
    ax = fig.add_axes((0.52, 0.58, 0.45, 0.4))
    profile = np.sum(gris_cutout*ok, axis=1)
    #ax.plot(profile, color='black')
    xx = np.arange(len(profile))
    ax.fill_between(xx-21, profile, profile*0, color='black', alpha=0.2)
    for i in range(5):
        profile = np.sum(models[i]*ok, axis=1)
        ax.plot(xx-21, profile, label=labels[i])
    
    ax.set_yticklabels([])
    #ax.set_xlim(6,36)
    ax.set_xlim(-15,15)
    #ax.set_xticks(np.arange(6,37,5))
    #ax.set_xticklabels(np.arange(6,37,5)-21)
    
    ax.set_xlabel(r'$y$ pixel')
    
    return fig
    
    
def compare():
    """
    Compare residuals of Nor's and my methods
    """
    import os
    import pyfits
    import numpy as np
    
    os.chdir('/Users/brammer/WFC3/STGrism/STGrismData')
    
    flt = pyfits.open('PREP_FLT/ibhm52c4q_flt.fits')
    flt_d = pyfits.open('PREP_FLT/ibhm52c2q_flt.fits')
    
    #gbb = pyfits.open('PREP_FLT/ibhm52c4q_flt_CONT.fits')
    gbb = pyfits.open('PREP_FLT/ibhm52c4q_flt_CONT_2x2.fits')
        
    nor_gauss140 = pyfits.open('Nor/1501109_0023102_149_G141.GaussF140W/FLT/ibhm52c4q_flt_2.CONT.fits')
    nor_gaussAll = pyfits.open('Nor/1501109_0023102_149_G141.GaussAll/FLT/ibhm52c4q_flt_2.CONT.fits')
    nor_seg140 = pyfits.open('Nor/1501109_0023102_149_G141.SegAll/FLT/ibhm52c4q_flt_2.CONT.fits')
    nor_segAll = pyfits.open('Nor/1501109_0023102_149_G141.SegF140W/FLT/ibhm52c4q_flt_2.CONT.fits')
    
    import threedhst.dq
    ds9 = threedhst.dq.myDS9()
    
    bad = (flt_d['DQ'].data & (4+4096)) > 0
    flt_d['SCI'].data[bad] = 0
    
    bad = (flt['DQ'].data & (4+4096)) > 0
    flt['SCI'].data[bad] = 0
    
    ds9.frame(1)
    ds9.view(flt[1].data)

    ds9.frame(2)
    ds9.view(flt[1].data-nor_gauss140[1].data)

    ds9.frame(3)
    ds9.view(flt[1].data-nor_gaussAll[1].data)

    ds9.frame(4)
    ds9.view(flt[1].data-nor_seg140[1].data)

    ds9.frame(5)
    ds9.view(flt[1].data-nor_segAll[1].data)

    ds9.frame(6)
    ds9.view(flt[1].data-gbb[0].data)
    
    ds9.set('pan to 442 398')
    ds9.set('zoom to 1.2')
    