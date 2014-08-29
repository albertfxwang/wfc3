"""
Wavelength calibration
"""
import numpy as np
import matplotlib.pyplot as plt

def check_vy22():
    """
    See if I can reproduce the Vy2-2 lines with the interlacing code, interlaced 1x1 with a single image
    """
    import mywfc3.wave
    
    os.chdir('/Users/brammer/WFC3/Calibration/Cycle20/13093_Wavelength/Reduce')
    
    asn = threedhst.utils.ASNFile('dummy_asn.fits')
    
    root=''
    images = {'F098M':'ic5v02a8q_flt.fits.gz', 'F105W':'ic5v02a9q_flt.fits.gz', 'G102':'ic5v02aaq_flt.fits.gz', 'F140W':'ic5v41awq_flt.fits.gz', 'F160W':'ic5v41axq_flt.fits.gz', 'G141':'ic5v41ayq_flt.fits.gz'}

    root='-off+x'
    images = {'F098M':'ic5v02afq_flt.fits.gz', 'F105W':'ic5v02agq_flt.fits.gz', 'G102':'ic5v02ahq_flt.fits.gz', 'F140W':'ic5v41b3q_flt.fits.gz', 'F160W':'ic5v41b4q_flt.fits.gz', 'G141':'ic5v41b5q_flt.fits.gz'}

    root='-off-x'
    images = {'F098M':'ic5v02aiq_flt.fits.gz', 'F105W':'ic5v02ajq_flt.fits.gz', 'G102':'ic5v02akq_flt.fits.gz', 'F140W':'ic5v41b6q_flt.fits.gz', 'F160W':'ic5v41b7q_flt.fits.gz', 'G141':'ic5v41b8q_flt.fits.gz'}

    blue = ['F098M', 'F105W', 'G102']
    
    flat_file = {'G102':os.getenv('iref')+'/uc72113oi_pfl.fits', #F105W
                 'G141':os.getenv('iref')+'/uc721143i_pfl.fits'} #F140W
    
    flat = {}
    for key in flat_file.keys():
        im = pyfits.open(flat_file[key])
        flat[key] = im[1].data[5:-5, 5:-5]
    
    for filter in images.keys():
        test = filter in blue
        band = 'blue'*test + 'red'*(not test)
        asn.product = 'Vy22%s-%s-%s' %(root, band, filter)
        asn.exposures = [images[filter].split('_flt')[0]]
        asn.write(asn.product+'_asn.fits')
        im = pyfits.open('../RAW/'+images[filter])
        if filter in flat.keys():
            im[1].data /= flat[key]
            sky = pyfits.open('/Users/brammer/3DHST/Spectra/Work/CONF/sky.G141.set002.fits ')[0].data
            #sky /= flat[key]
            ratio = im[1].data/sky
            #a = np.median(ratio)
            #a = np.median(ratio[ratio < a*1.5])
            yh, xh = np.histogram(ratio.flatten(), range=(0,10), bins=1000)
            a = xh[1:][np.argmax(yh)]
            bg = a
            im[1].data -= a*sky
        else:
            bg = np.median(im[1].data)
            im[1].data -= bg
        #
        print 'Background: %s %.4f' %(images[filter], bg)
        #
        #if not os.path.exists(images[filter][:-3]):
        im.writeto(images[filter].split('.gz')[0], clobber=True)
    
    files=glob.glob('Vy22%s-[br]*asn.fits' %(root))
    for file in files:
        unicorn.reduce.interlace_combine(file.split('_asn')[0], growx=1, growy=1, pad=60, NGROW=100, view=False)
    
    #### determine shifts to make spectra smooth at the edges
    # shifts = {'Vy22-red-G141':(0,1), 'Vy22-blue-G102':(0,1)}
    # shifts = {'Vy22-off+x-red-G141':(0,1), 'Vy22-off+x-blue-G102':(0,1)}
    # shifts = {'Vy22-off-x-red-G141':(0,1), 'Vy22-off-x-blue-G102':(0,1)}
    # for root in shifts.keys():
    #     im = pyfits.open(root+'_inter.fits', mode='update')
    #     for ext in [1,2]:
    #         for axis in [0,1]:
    #             im[ext].data = np.roll(im[ext].data, shifts[root][axis], axis=axis)
    #     #
    #     im.flush()
        
    fig = unicorn.plotting.plot_init(xs=10, aspect=0.5, left=0.1, bottom=0.1, wspace=0, hspace=0)
    
    ### Run twice with old and new configuration files
    new, sub = False, 211
    new, sub = True, 212
    
    ax = fig.add_subplot(sub)

    for root, alpha in zip(['', '-off-x', '-off+x'], [0.4,0.4,0.4]):
    #for root, alpha in zip(['', '-off-x', '-off+x'][:1], [0.4,0.4,0.4][:1]):
        sp = mywfc3.wave.get_vy22(root='Vy22%s-blue' %(root), new=new)
        ax.plot(sp.oned.lam, sp.oned.flux, color='blue', linewidth=2, alpha=alpha)
        sp = mywfc3.wave.get_vy22(root='Vy22%s-red' %(root), new=new)
        ax.plot(sp.oned.lam, sp.oned.flux, color='red', linewidth=2, alpha=alpha)

    ax.semilogy()
    
    PNe_lines = [9071.403457, 9534.921052, 10049.850283, 10833.000000, 12821.000000, 16112.000000, 16412.000000]
    ## Paper
    #PNe_lines = [11621, 11665, 11892, 11970, 12529, 12817, 15335, 15549, 15693, 15875, 16102, 16401, 16801]
    for line in PNe_lines:
        ax.plot(np.ones(2)*line, [100,1.e5], color='black', alpha=0.5, linewidth=2)
    
    for ax in fig.axes:
        ax.set_xlim(7400,1.68e4)
        ax.set_ylim(300, 4.e4)
    
    ax.set_xlabel(r'$\lambda$')
    
    unicorn.plotting.savefig(fig, 'Vy22_center_Fixed_v2.pdf')
    
    ##### Full model
    root='Vy22-red'
    
    if 'blue' in root:
        direct='F105W'
        grism='G102'
    #
    if 'red' in root:
        direct='F140W'
        grism='G141'
    
    model = unicorn.reduce.GrismModel(root, direct=direct, grism=grism, growx=1, growy=1)
    model.compute_full_model(BEAMS=['A', 'B', 'C', 'D', 'E'], view=None, MAG_LIMIT=18.0, save_pickle=True, refine=False, model_slope=-0.5)

    model = unicorn.reduce.GrismModel(root, direct=direct, grism=grism, growx=1, growy=1)
    model.compute_full_model(BEAMS=['B'], view=None, MAG_LIMIT=20.0, save_pickle=True, refine=False, model_slope=-0.5)
    
    sp = unicorn.reduce.Interlace2D('Vy22-red_00602.2D.fits')
    plt.plot(sp.oned.lam, sp.oned.flux)
    yi = np.interp(1.4e4, sp.oned.lam, sp.oned.flux)
    plt.plot(sp.oned.lam, yi*(sp.oned.lam/1.4e4)**beta)
    
    im = pyfits.open('%s-G141_inter.fits' %(root))
    

def brightest(root='Vy22-red', N=10):
    new=False
    
    import numpy as np
    import unicorn
    
    root='Vy22-red'
    N=10
    
    if 'blue' in root:
        direct='F105W'
        grism='G102'
    #
    if 'red' in root:
        direct='F140W'
        grism='G141'
    #
    NGROW=135
    files=glob.glob('%s-*asn.fits' %(root))
    for file in files:
        unicorn.reduce.interlace_combine(file.split('_asn')[0], growx=1, growy=1, pad=60, NGROW=NGROW, view=False)
    #
    os.system('rm *seg.fits *cat')
    
    unicorn.reduce.set_grism_config(grism=grism, chip=1, use_new_config=new, force=False)
    model = unicorn.reduce.GrismModel(root, direct=direct, grism=grism, growx=1, growy=1)
    
    objects = model.objects[np.argsort(model.cat.mag)][:N]
    for object in objects:
        model.twod_spectrum(object, miny=-30, refine=False)
        
    objects = model.objects[np.argsort(model.cat.mag)][:N]
    for object in objects:
        sp = unicorn.reduce.Interlace2D('%s_%05d.2D.fits' %(root, object))
        wave, flux = sp.optimal_extract(sp.im['SCI'].data-sp.im['CONTAM'].data)
        sens = sp.im['SENS'].data    
        yi = np.interp(1.5e4, wave, flux/sens)
        plt.plot(wave, flux/sens/yi, color='blue', alpha=0.5, linewidth=3)
        
def get_vy22(root='Vy22-blue', new=False):
    import numpy as np
    import unicorn
    
    if 'blue' in root:
        direct='F105W'
        grism='G102'
    #
    if 'red' in root:
        direct='F140W'
        grism='G141'
    
    unicorn.reduce.set_grism_config(grism=grism, chip=1, use_new_config=new, force=False)
    model = unicorn.reduce.GrismModel(root, direct=direct, grism=grism, growx=1, growy=1)
    
    r0, d0 =  291.0926, 9.8989242
    bright = model.cat.mag < 12.2
    dr = np.sqrt(((model.ra_wcs[bright]-r0)*np.cos(d0/180*np.pi))**2+(model.dec_wcs[bright]-d0)**2)*3600
    id = model.objects[bright][dr == dr.min()][0]
    
    print root, id, dr.min()
    model.twod_spectrum(id, miny=-20, refine=False)
    sp = unicorn.reduce.Interlace2D('%s_%05d.2D.fits' %(root, id))
    return sp
      
def find_zeroth():
    root='Vy22-red'
    if 'blue' in root:
        direct='F105W'
        grism='G102'
    #
    if 'red' in root:
        direct='F140W'
        grism='G141'
    
    unicorn.reduce.set_grism_config(grism=grism, chip=1, use_new_config=False, force=True)
    model = unicorn.reduce.GrismModel(root, direct=direct, grism=grism, growx=1, growy=1)
    
    im = pyfits.open(root+'_inter_model.fits')
    resid = model.gris[1].data-im[0].data*0
    
    so = np.argsort(model.cat.mag)
    
    r = np.sqrt((model.cat.x_pix-567)**2 + (model.cat.y_pix-567)**2)
    so = np.argsort(r)[::-1]
    
    pad = 60
    
    N = 4
    yc, xc = np.indices((2*N, 2*N))
    
    xin = []
    yin = []
    yout = []
    
    for i in so[model.cat.mag[so] < 17]:
        x0, y0 = mywfc3.wave.predict_zeroth(model.cat.x_pix[i]-pad, model.cat.y_pix[i]-pad)
        x0i, y0i = int(x0)+pad, int(y0)+pad
        if (x0i < pad) | (y0i < pad):
            continue
        #
        sub = resid[y0i-N:y0i+N, x0i-N:x0i+N]
        ok = sub > 0
        mask = sub == 0
        med = nd.median_filter(sub, 3)
        sub[mask] = med[mask]
        yf0 = np.sum((yc*sub))/np.sum(sub) 
        yf = yf0 - N + y0i - pad + 1
        xf0 = np.sum((xc*sub))/np.sum(sub) 
        xf = xf0 - N + x0i - pad + 1
        xin.append([model.cat.x_pix[i]-pad, model.cat.y_pix[i]-pad])
        yin.append([x0, y0])
        yout.append([xf, yf])
        ds9.view(sub)
        
    xin = np.array(xin).T
    yin = np.array(yin).T
    yout = np.array(yout).T

def ic5117():
    """
    Check grism overlap with IC-5117 observations (13582)
    """
   import threedhst
   import numpy as np
   import glob
   import os
   import pyfits
   os.chdir('/Users/brammer/WFC3/Calibration/Cycle21/13582_Wavelength/PREP_FLT')
    
    asn = threedhst.utils.ASNFile('dummy_asn.fits')
    files=glob.glob('*flt.fits')
    for file in files:
        asn.exposures = [file.split('_flt')[0]]
        im = pyfits.open(file, mode='update')
        filt=im[0].header['FILTER']
        iref = os.getenv('iref')
        if filt in ['F098M','F127M','F105W','G102']:
            gr = 'b'
            flat = pyfits.open(os.path.join(iref, 'uc72113oi_pfl.fits'))
        else:
            gr = 'r'
            flat = pyfits.open(os.path.join(iref, 'uc721143i_pfl.fits'))
        #
        if 'FLAT' not in im[0].header.keys():
            im[0].header['FLAT'] = os.path.basename(flat.filename())
            im[1].data /= flat[1].data[5:-5,5:-5]
            im[1].data -= np.median(im[1].data)
            im.flush()
        #
        asn.product = 'IC5117%s-%s' %(gr, filt)
        asn.write('%s_asn.fits' %(asn.product))
        
    
    unicorn.reduce.interlace_combine('IC5117b-F105W', view=False, growx=1, growy=1, auto_offsets=True)
    unicorn.reduce.interlace_combine('IC5117b-F098M', view=False, growx=1, growy=1, auto_offsets=True)
    unicorn.reduce.interlace_combine('IC5117b-G102', view=False, growx=1, growy=1, auto_offsets=True)

    unicorn.reduce.interlace_combine('IC5117r-F140W', view=False, growx=1, growy=1, auto_offsets=True)
    unicorn.reduce.interlace_combine('IC5117r-G141', view=False, growx=1, growy=1, auto_offsets=True)
    
    model = unicorn.reduce.GrismModel('IC5117r', grism='G141', direct='F140W', growx=1, growy=1)
    model.twod_spectrum(114, miny=-80, refine=False, extract_1d=True)
    model.show_2d(savePNG=True)
    
    os.system('ln -s IC5117r_inter_seg.fits IC5117b_inter_seg.fits')
    os.system('ln -s IC5117r_inter.reg IC5117b_inter.reg')
    os.system('ln -s IC5117r_inter.cat IC5117b_inter.cat')

    model = unicorn.reduce.GrismModel('IC5117b', grism='G102', direct='F105W', growx=1, growy=1)
    model = unicorn.reduce.GrismModel('IC5117b', grism='G102', direct='F098M', growx=1, growy=1)
    model.twod_spectrum(114, miny=-80, refine=False, extract_1d=True)
    model.show_2d(savePNG=True)
    
    #
    gris = unicorn.interlace_fit.GrismSpectrumFit('IC5117b_00114', skip_photometric=True)
    rx, ry = np.loadtxt('rudy_spec.dat', unpack=True)
    gris.twod.compute_model(rx+7, ry/1.e-17/gris.twod.total_flux) #, norm=False)

    gris = unicorn.interlace_fit.GrismSpectrumFit('IC5117r_00114', skip_photometric=True)
    rx, ry = np.loadtxt('rudy_spec.dat', unpack=True)
    gris.twod.compute_model(rx-20, ry/1.e-17/gris.twod.total_flux) #, norm=False)
    
    twod = gris.twod    
    w, f = twod.optimal_extract(twod.im['SCI'].data)
    plt.plot(w, f/twod.im['SENS'].data, color='black', alpha=0.8, linewidth=2)
    w, f = twod.optimal_extract(twod.model)
    plt.plot(w, f/twod.im['SENS'].data, color='red', alpha=0.8, linewidth=2)
    #plt.plot(nd.shift(w, 1), (f/twod.im['SENS'].data), color='orange', alpha=0.8, linewidth=2)
    
    plt.xlim(8000,1.6e4)
    plt.semilogy()
    plt.ylim(500, 3.5e4)
    
    t = table.read('rudy_table1.dat.txt', format='ascii.cds')
    
    cont = S.FlatSpectrum(6e-11*1.e-4*1.e-8, fluxunits='photlam')
    
    templ = np.ones((len(t), len(w)))
    use = np.zeros(len(t))
    
    for i in range(len(t)):
        line = t[i]
        if line['f_Ratio'] == 'Atmosphere':
            continue
        #
        if line['Ratio'] < 0.02:
            continue
        #
        use[i] = 1
        spec = cont + S.GaussianSource(line['Ratio']*4.95e-12, line['Wave']*1., line['Wave']*800./3.e5, fluxunits='photlam')
        gris.twod.compute_model(spec.wave, spec.flux/1.e-17/gris.twod.total_flux) #, norm=False)
        w, f = twod.optimal_extract(twod.model)
        templ[i,:] = f
        
    ok = use > 0
    templ = templ[ok,:]
    cont = S.FlatSpectrum(6e-11*1.e-4, fluxunits='photlam')
    gris.twod.compute_model(cont.wave, cont.flux/1.e-17/gris.twod.total_flux) #, norm=False)
    w, f = twod.optimal_extract(twod.model)
    templ[0,:] = f
    
    wobs, fobs = twod.optimal_extract(twod.im['SCI'].data)
    var = 1./twod.im['SENS'].data
    var = (var/var.min()*0.01*fobs.max())**2
    var = (var/var.min()*0.001*fobs.max())**2 ### for blue
    
    amatrix = unicorn.utils_c.prepare_nmf_amatrix(var, templ)
    coeff = unicorn.utils_c.run_nmf(fobs, var, templ, amatrix, verbose=True, toler=1.e-10)
    plt.plot(wobs, fobs/twod.im['SENS'].data, color='black', alpha=0.8, linewidth=2)
    sum = np.dot(coeff, templ)
    plt.plot(wobs, sum/twod.im['SENS'].data, color='red', alpha=0.8, linewidth=2)
    sum2 = np.dot(coeff*0+1, templ)
    plt.plot(wobs, sum2/twod.im['SENS'].data, color='green', alpha=0.8, linewidth=2)
    
    #sum = np.dot(coeff*0+1, templ)
    
def rudy_ic5117():
    """
    Make a simple spectrum from the line table from Rudy et al.
    """
    from astropy.table import Table as table
    import pysynphot as S
    t = table.read('rudy_table1.dat.txt', format='ascii.cds')
    
    spec = S.FlatSpectrum(6e-11*1.e-4, fluxunits='photlam')
    #spec.convert('photlam')
    for line in t:
        if line['f_Ratio'] == 'Atmosphere':
            continue
        #
        #spec += S.GaussianSource(line['Ratio']*4.95e-12, line['Wave']*1., line['Wave']*30./3.e5, fluxunits='photlam')
        if line['Wave'] == 10830:
            f = 1.35
        else:
            f = 1.
        #
        spec += S.GaussianSource(line['Ratio']*4.95e-12*f, line['Wave']*1., line['Wave']*800./3.e5, fluxunits='photlam')
    
    ok = (spec.wave > 7000) & (spec.wave < 1.7e4)
    np.savetxt('rudy_spec.dat', np.array([spec.wave[ok], spec.flux[ok]]).T, fmt='%.5e')
    
def predict_zeroth(xc, yc, new=False, grism='G141', force=False):
    import unicorn
    
    unicorn.reduce.set_grism_config(grism=grism, chip=1, use_new_config=new, force=force)
    
    beam = 'B'
    conf = unicorn.reduce.conf
    
    growx, growy = 1,1
    
    bigX = xc
    bigY = yc
    
    dydx_order = np.int(conf['DYDX_ORDER_'+beam])
    dydx_0 = unicorn.reduce.field_dependent(bigX, bigY, conf['DYDX_'+beam+'_0']) * growy
    dydx_1 = unicorn.reduce.field_dependent(bigX, bigY, conf['DYDX_'+beam+'_1']) * growy / growx#/ grow_factor        
                
    xoff_beam = unicorn.reduce.field_dependent(bigX, bigY, conf['XOFF_'+beam]) * growx
    yoff_beam = unicorn.reduce.field_dependent(bigX, bigY, conf['YOFF_'+beam]) * growy
    
    disp_order = np.int(conf['DISP_ORDER_'+beam])
    dldp_0 = unicorn.reduce.field_dependent(bigX, bigY, conf['DLDP_'+beam+'_0'])
    dldp_1 = unicorn.reduce.field_dependent(bigX, bigY, conf['DLDP_'+beam+'_1']) / growx
    
    #lam = dldp_0 + dldp_1*(xarr-xoff_beam)
    
    x0 = (1.4e4 - dldp_0)/dldp_1 + xoff_beam
    return xc + x0, yc + yoff_beam
    
def offset_position():
    """
    Check grism extractions for 13580 program that put Vy2-2 just
    off of the left edge of the frame to test extrapolation of wavelength
    calibration.
    """
    
    import astropy.io.fits as pyfits
    from astropy.table import Table as table
    
    import drizzlepac
    from drizzlepac import tweakreg, tweakback
    import stwcs
    
    import unicorn
    
    unicorn.candels.make_asn_files(uniquename=True)
    
    info = table.read('files.info', format='ascii.commented_header')
    
    for filter in ['F098M', 'F105W']:
        filter_files = list(info['FILE'][info['FILTER'] == filter])
        #
        files = glob.glob('VY2-2*%s_asn.fits' %(filter))
        for file in files:
            prep.prep_direct_grism_pair(direct_asn=file, grism_asn=False, radec='2mass.radec', scattered_light=False, skip_direct=False)  
        #
        driz_images = glob.glob('VY2-2*%s_drz_sci.fits' %(filter))
        tweakreg.TweakReg(driz_images, refimage=driz_images[0], updatehdr=True, updatewcs=True, catfile=None, xcol=2, ycol=3, xyunits='pixels', refcat=None, refxcol=1, refycol=2, refxyunits='degrees', shiftfile=True, outshifts='%s_shifts.txt' %(filter), outwcs='%s_wcs.fits' %(filter), searchrad=5, tolerance=12, wcsname='TWEAK', interactive=False, residplot='No plot', see2dplot=False, clean=True, headerlet=True, clobber=True)
        tweakback.tweakback(driz_images[1])
        #
        drizzlepac.astrodrizzle.AstroDrizzle(filter_files, output='VY22-%s' %(filter),  clean=True, skysub=False, final_scale=None, final_pixfrac=1, context=False, final_bits=576, preserve=False, driz_cr_snr='5.0 4.0', driz_cr_scale = '2.5 0.7')
        drizzlepac.astrodrizzle.AstroDrizzle(filter_files, output='VY22-%s' %(filter), clean=True, context=False, preserve=False, skysub=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True)
        
    ### Put WCS from direct F105W images into G102 at same POS-TARG
    info = table.read('files.info', format='ascii.commented_header')
    
    idx = np.arange(len(info))[info['FILTER'] == 'F105W']
    asn = threedhst.utils.ASNFile('../RAW/ibhj01030_asn.fits')
    
    for i in idx:
        direct = info['FILE'][i]
        dx, dy = info['POSTARG1'][i], info['POSTARG2'][i]
        ix_gris = (info['POSTARG1'] == dx) & (info['POSTARG2'] == dy) & (info['FILTER'] == 'G102')
        grism = info['FILE'][ix_gris][0]
        sign = {True:'+', False:'-'}
        #
        asn.product = 'VY22%s%02d%s%02d-F105W' %(sign[dx > 0], np.abs(dx), sign[dy > 0], np.abs(dy))
        asn.exposures = [direct.split('_flt')[0]]
        asn.write(asn.product + '_asn.fits')
        #
        asn.product = 'VY22%s%02d%s%02d-G102' %(sign[dx > 0], np.abs(dx), sign[dy > 0], np.abs(dy))
        asn.exposures = [grism.split('_flt')[0]]
        asn.write(asn.product + '_asn.fits')
        #### update WCS header
        imd = pyfits.open(direct)
        img = pyfits.open(grism)
        sci_ext=1
        direct_WCS = stwcs.wcsutil.HSTWCS(imd, ext=sci_ext)
        drizzlepac.updatehdr.update_wcs(grism, sci_ext, direct_WCS, verbose=True)    
        
    #### Make reference catalog
    root = 'VY22-F105W'
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION, BACKGROUND'
    se.options['CHECKIMAGE_NAME'] = '%s_drz_seg.fits, %s_drz_bkg.fits' %(root, root)
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = '%s_drz_wht.fits' %(root)
    se.options['WEIGHT_GAIN'] = 'Y'
    se.options['GAIN'] = '0'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_THRESH']    = '2.' 
    se.options['ANALYSIS_THRESH']  = '2.' 
    se.options['DETECT_MINAREA'] = '10' 
    se.options['MASK_TYPE'] = 'NONE'
    se.options['DEBLEND_NTHRESH']    = '64' 
    se.options['DEBLEND_MINCONT']  = '0.1' 
    se.options['SEEING_FWHM'] = '0.12'
    
    se.options['BACK_TYPE'] = 'MANUAL'
    se.options['BACKPHOTO_TYPE'] = 'LOCAL'
    
    se.options['MAG_ZEROPOINT'] = '%.2f' %(unicorn.reduce.ZPs['F105W'])
    se.options['CATALOG_TYPE']    = 'ASCII_HEAD'
    se.options['CATALOG_NAME']    = '%s_drz_sci.cat' %(root)
    status = se.sextractImage('%s_drz_sci.fits[0]' %(root))
    threedhst.sex.sexcatRegions('%s_drz_sci.cat' %(root), '%s_drz_sci.reg' %(root), format=1)
    
    #### Make interlaced images
    files = glob.glob('VY22[+-]??[+-]??-F105W_asn.fits')
    for file in files:
        unicorn.reduce.interlace_combine(file.split('_asn')[0], growx=1, growy=1, NGROW=50, pad=60, view=False)
        unicorn.reduce.interlace_combine(file.split('_asn')[0].replace('F105W', 'G102'), growx=1, growy=1, NGROW=50, pad=60, view=False)
        red = unicorn.reduce
        red.adriz_blot_from_reference(pointing=file.split('_asn')[0], pad=60, NGROW=50, growx=1, growy=1, auto_offsets=False, ref_exp=0, ref_image='VY22-F105W_drz_sci.fits', ref_ext=0, ref_filter='F105W', seg_image='VY22-F105W_drz_seg.fits', cat_file='VY22-F105W_drz_sci.cat')
    
    ### extract spectra
    id = 798
    files = glob.glob('VY22[+-]??[+-]??-F105W_asn.fits')
    for file in files:
        model = unicorn.reduce.GrismModel(root=file.split('-F10')[0], direct='F105W', grism='G102', growx=1, growy=1, grow_factor=1)
        model.twod_spectrum(id, miny=-30, refine=False, CONTAMINATING_MAGLIMIT=0)
        
    files = glob.glob('*2D.fits')
    yi, xi = np.indices((30,30))
    xs, ys = np.zeros(len(files)), np.zeros(len(files))
    xpix = xs*0
    for i, file in enumerate(files):
        twod = unicorn.reduce.Interlace2D(file)
        xs[i] = np.sum(xi*twod.im['DSCI'].data/twod.im['DSCI'].data.sum())
        ys[i] = np.sum(xi*twod.im['DSCI'].data/twod.im['DSCI'].data.sum())
        xpix[i] = twod.im[0].header['X_PIX']
        
    xs -= np.median(xs) #+ 0.5
    ys -= np.median(ys)
    
    #xs -= 0.5
    
    fig = plt.figure(figsize=[16,4])
    fig.subplots_adjust(left=0.04, right=0.98, top=0.92)
    for i, file in enumerate(files):
        twod = unicorn.reduce.Interlace2D(file)
        w, f = twod.optimal_extract(twod.im['SCI'].data)
        c = {True: 'red', False: 'blue'}
        #plt.plot(w-np.diff(w)[0]*xs[i], f/twod.im['SENS'].data, alpha=0.5, marker='o', ms=2, label='%s, %s' %(file[4:7], file[7:10])) # , color=c['2-2' in file]
        ff = f*0.
        for k in range(ff.shape[0]):
            y0 = int(np.round(twod.im['YTRACE'].data[k]))
            ff[k] = np.sum(twod.im['SCI'].data[y0-4:y0+4, k])
        #
        plt.plot(w-np.diff(w)[0]*xs[i], ff/twod.im['SENS'].data, alpha=0.5, marker='o', ms=2, label='%s, %s' %(file[4:7], file[7:10])) # , color=c['2-2' in file]
        #plt.plot(twod.oned.data['wave']-np.diff(w)[0]*xs[i], twod.oned.data['flux']/twod.oned.data['sensitivity'], alpha=0.5, marker='o', ms=2, label='%s, %s' %(file[4:7], file[7:10])) # , color=c['2-2' in file]
        #
        print file, np.diff(w)[0]
        #ds9.frame(i+1)
        #ds9.view(twod.im['DSCI'].data)
    
    PNe_lines = [9071.403457, 9534.921052, 10049.850283, 10833.000000, 12821.000000, 16112.000000, 16412.000000]
    for line in PNe_lines:
        plt.plot([line, line], [0.1,1.e5], color='black', linewidth=3, alpha=0.2, zorder=-5)
        
    #plt.plot(w-np.diff(w)[0]*xs[i]-np.diff(w)[0], f/twod.im['SENS'].data, alpha=0.5, color='green', marker='o', ms=2)
    plt.legend(loc='upper right', prop={'size':9}, title='POS-TARG')
    plt.title('VY2-2, G102, 13580')
    plt.xlim(8500, 11500)
    plt.ylim(700,14000)
    plt.ylim(600,64000)
    plt.semilogy()
    plt.xlabel(r'$\lambda$')
    plt.savefig('vy22-edge_v2.pdf') #, dpi=100)
    
    plt.close()
    