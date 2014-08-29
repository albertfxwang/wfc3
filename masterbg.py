"""
Make master grism backgrounds from masked FLT images
"""

def make_ima():
    import glob
    from subprocess import Popen,PIPE
    import os

    import wfc3tools
      
    #os.chdir("/Volumes/3DHST_Gabe/IMA")
    os.chdir('/Volumes/3DHST_Gabe/MasterBackgroundG102/IMA')
    grism = 'G102'
    grism = 'G141'
    
    files=glob.glob('*raw.fits')
    for file in files: 
        stdout, stderr = Popen('gethead -x 0 %s FILTER' %(file), shell=True, stdout=PIPE).communicate()    
        #print file
        if stdout.split()[0] != grism:
            print stdout
            root=file.split('q_')[0]
            rems = glob.glob('%s*' %(root))
            for rem in rems:
                print rem
                os.remove(rem)
            #
            continue
        #
        if os.path.exists(file.replace('raw','ima')):
            continue
        #
        try:
            wfc3tools.calwf3.calwf3(file)
            os.remove(file.replace('raw','flt'))
        except:
            pass
        
    #### Copy flt.seg.fits files
    GRISM = 'G141'
    files=glob.glob('ic[3d]*%s_orbit.dat' %(GRISM))
    for file in files:
        mask_file = file.replace('j_%s_orbit.dat' %(GRISM), 'q_flt.seg.fits')
        if not os.path.exists(mask_file):
            PATH = '/user/brammer/WFC3_Backgrounds/WISPS/%s/' %(GRISM)
            PATH = '/user/brammer/WFC3_Backgrounds/WISPS/NewG141/%s' %(GRISM)
            print mask_file
            print os.path.exists('%s/%s' %(PATH, mask_file))
            os.system('rsync -avz %s/%s ./' %(PATH, mask_file))
            
def mask_all():
    os.chdir("/Users/brammer/WFC3/Backgrounds/MaskedMaster/G141")
    files=glob.glob('*model.fits')
    for file in files:
        root=file.split('_inter')[0]
        mywfc3.masterbg.mask_visit(root=root, threshold=0.005, grow_mask=8)
    
def make_cleaned():
    """
    Make images of masked individual reads
    """   
    import numpy.ma as ma
    import pyfits
    import unicorn
    from threedhst import catIO
    import os
    import glob
    import scipy.ndimage as nd
    
    #flat_f140 = pyfits.open('%s/flat_3DHST_F140W_t2_v0.1.fits' %(os.getenv('iref')))[1].data#[5:-5,5:-5]
    flat_f140 = pyfits.open('%s/uc721143i_pfl.fits' %(os.getenv('iref')))[1].data#[5:-5,5:-5]
    flat_f105 = pyfits.open('%s/uc72113oi_pfl.fits' %(os.getenv('iref')))[1].data#[5:-5,5:-5]
    flat_g141 = pyfits.open('%s/u4m1335mi_pfl.fits' %(os.getenv('iref')))[1].data#[5:-5,5:-5]
    flat_g102 = pyfits.open('%s/u4m1335li_pfl.fits' %(os.getenv('iref')))[1].data
        
    #flat = pyfits.open('%s/flat_3DHST_F140W_t2_v0.1.fits' %(os.getenv('iref')))[1].data #[5:-5,5:-5]
    flat = pyfits.open('%s/uc721143i_pfl.fits' %(os.getenv('iref')))[1].data
    flat_g141 = pyfits.open('%s/u4m1335mi_pfl.fits' %(os.getenv('iref')))[1].data
    flat /= flat_g141

    PATH = '/Users/brammer/WFC3/Backgrounds/MaskedMaster/Go_G141'
    zodi = pyfits.open('%s/zodi_G141_clean.fits' %(PATH))[0].data
    zodi_wht = pyfits.open('%s/zodi_G141_clean_wht.fits' %(PATH))[0].data
    zodi_dq = pyfits.open('%s/zodi_G141_clean_dq.fits' %(PATH))[0].data
    zodi_flat = pyfits.open('%s/zodi_G141_flat.fits' %(PATH))
    flat = zodi_flat[1].data/flat_g141

    ### G102
    os.chdir("/Volumes/3DHST_Gabe/MasterBackgroundG102")
    mask_ext = '_flt.seg'
    files = glob.glob('IMA/*%s.fits' %(mask_ext))
    GRISM = 'G102'
    
    flat = flat_f105/flat_g102
    PATH = '/Users/brammer/WFC3/Backgrounds/MaskedMaster/G102'
    PATH = '/Volumes/3DHST_Gabe/MasterBackgroundG102'
    zodi = pyfits.open('%s/zodi_G102_clean.fits' %(PATH))[0].data
    zodi_dq = zodi*0
    
    ### G141
    os.chdir("/Volumes/3DHST_Gabe/MasterBackground")
    mask_ext = '_flt_mask'
    #mask_ext = '_flt.seg' ### WISPS
    files = glob.glob('IMA/*%s.fits' %(mask_ext))
    GRISM = 'G141'
    
    flat = flat_f140/flat_g141
    #PATH = '/Users/brammer/WFC3/Backgrounds/MaskedMaster/Go_G141'
    PATH = '/Volumes/3DHST_Gabe/MasterBackground'
    zodi = pyfits.open('%s/zodi_G141_clean.fits' %(PATH))[0].data
    zodi_dq = zodi*0
    
    SKIP = True
    
    N = 0
    
    for file in files:
        if not os.path.exists(file.replace(mask_ext, '_ima')):
            continue
        if not os.path.exists(file.split('q_')[0]+'j_%s_orbit.dat' %(GRISM)):
            continue
        #
        zodi_level, zodi_all = mywfc3.masterbg.visit_minimum(base=file, PATH='./IMA/')
        ### G141
        root = os.path.basename(file.split(mask_ext)[0])
        bg = catIO.Readfile(file.split('q_')[0] + 'j_%s_orbit.dat' %(GRISM))
        #
        level, ok_reads = 'low', np.arange(bg.N)[bg.shadow > 0] + 1
        ###level, ok_reads = 'high', np.arange(bg.N)[(bg.bg/zodi_level > 2.0) & (bg.shadow == 0)] + 1
        level, ok_reads = 'high', np.arange(bg.N)[(bg.bg-zodi_level > 0.7) & (bg.shadow == 0)] + 1
        N+= len(ok_reads)
        print N
        if ok_reads.sum() == 0:
            continue
        #
        if (len(glob.glob('%s*%s*' %(root, level))) > 0) & SKIP:
            print root
            continue
        else:
            print 'Run: %s' %(root)
        #
        ima = pyfits.open('IMA/'+root+'_ima.fits')
        cube, dq, time, NSAMP = unicorn.prepare.split_multiaccum(ima)
        diff, dt = np.diff(cube, axis=0), np.diff(time)
        mask = pyfits.open(file)[0].data > 0
        #
        NGROW = 12
        mask = nd.maximum_filter(mask, size=NGROW)
        #
        NSAMP = ima[0].header['NSAMP']
        for read in ok_reads:
            head = ima['sci',NSAMP-(read+1)].header.copy()
            #im_read = (diff[read]/dt[read]/flat)[5:-5,5:-5]
            #scl = np.median(im_read[~mask])
            #im_read /= scl
            #im_read[mask | ((dq[read+1, :,:][5:-5,5:-5] & 32) == 32)] = -100
            ## subtract zodi BG for high-bg case
            im_read = (diff[read]/dt[read]/flat)[5:-5,5:-5]
            if level == 'high':
                im_read -= zodi*zodi_level
            #
            scl = np.median(im_read[~mask])
            im_read /= scl
            im_read[mask | ((dq[read+1, :,:][5:-5,5:-5] & (4+16+32)) > 0) | (zodi_dq > 0)] = -100
            head.update('BGSCL', scl)
            head.update('BGOBS', bg.bg[read-1])
            head.update('BGZODI', zodi_level)
            head.update('BGREAD', bg.read[read-1])
            head.update('BGLIMB', bg.limbang[read-1])
            out = '%s_%s_%02d.fits' %(root, level, read)
            #out = '%s_low_%02d.fits' %(root, read)
            pyfits.writeto(out, data=im_read, header=head, clobber=True)
            mx = ma.masked_array(im_read, mask=(im_read < -90))
            med = ma.median(mx, axis=0)
            #plt.plot(med, alpha=0.1, color='blue')
            pyfits.writeto(out.replace('.fits','_med.fits'), data=med.data, header=head, clobber=True)
            #
            print out
            
    ### Show 1D
    import os
    import glob
    import numpy.ma as ma
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors
    import pyfits
    
    os.chdir("/Volumes/3DHST_Gabe/MasterBackground")
    #os.chdir("/Volumes/3DHST_Gabe/MasterBackgroundG102")
    
    #os.chdir("/Users/brammer/WFC3/Backgrounds/MaskedMaster/Go_G141")
    level = 'low'
    level = 'high'
    files, color = glob.glob('i*%s*med.fits' %(level)), 'blue'
        
    N = len(files)
    xarr = np.arange(1014)
    #cm = matplotlib.colors.Normalize(vmin=0.5, vmax=4)
    yarr = np.zeros((len(files), 1014))
    xglint = np.zeros((len(files), 3))
    yglint = np.zeros(len(files))
    for i, file in enumerate(files[0:N]):
        medx = pyfits.open(files[i])
        print '%s  %.3f' %(file, np.std(medx[0].data))
        yarr[i,:] = medx[0].data*1.
        yarr[i,:] /= np.median(yarr[i,-350:-100])
        glint = np.median(yarr[i,200:400])/np.median(yarr[i,750:900])
        #xglint[i,:] = np.array([medx[0].header['BGLIMB'], medx[0].header['BGSCL'], medx[0].header['BGZODI']])
        xglint[i,:] = np.array([medx[0].header['BGOBS'], medx[0].header['BGSCL'], medx[0].header['BGZODI']])
        yglint[i] = glint
        medx.close()
        #print file
    
    skip = np.maximum(1, N/400)    
    for i in range(len(files))[0:N:skip]: #[::len(files)/200]:
        glint = np.median(yarr[i,200:400])/np.median(yarr[i,750:900])
        if glint > 1.05:
            continue
            plt.plot(xarr, yarr[i,:], alpha=0.1, color='red')
            print 'glint: %d, limbang: %f' %(i, xglint[i])
        else:
            plt.plot(xarr, yarr[i,:], alpha=0.01+0.1*(N < 100), color=color)
    
    mx = ma.masked_array(yarr)
    s = np.std(yarr, axis=0)
    m = np.median(yarr, axis=0)
    plt.plot(np.median(yarr, axis=0), color=color, linewidth=2, alpha=0.8)
    plt.plot(np.median(yarr, axis=0)+s, color='black', alpha=0.5)
    plt.plot(np.median(yarr, axis=0)-s, color='black', alpha=0.5)
    plt.plot(np.median(yarr, axis=0)+s, color=color, alpha=0.5)
    plt.plot(np.median(yarr, axis=0)-s, color=color, alpha=0.5)
    
    #### Color lines by some variable, say LimbAngle
    yshow = yarr*1
    for i in range(yarr.shape[0]):
        if (yglint[i] > 1.05) | (xglint[i,2] > 0.9):
            yshow[i,:] = xarr*0.+1
            
    sl = slice(0, yarr.shape[0], 1)
    segments = [list(zip(xarr[::4], yshow[i][::4])) for i in np.arange(yshow.shape[0])[sl]]
    import matplotlib.collections
    lc = matplotlib.collections.LineCollection(segments) #, norm=plt.Normalize(0.8,3)) #, cmap=plt.get_cmap('jet'))
    lc.set_array(xglint[sl, 0])
    lc.set_alpha(0.2)
    plt.plot(xarr, np.median(yarr, axis=0), color='black', zorder=10, linewidth=3, alpha=0.6)
    plt.gca().add_collection(lc)
    cb = plt.colorbar(lc)
    
    #### plot the column-averaged curves in a single image
    mmid = np.median(yarr[:,200:800], axis=1)
    mlow = np.median(yarr[:,10:84], axis=1)
    bad = (np.abs(mmid-np.median(mmid)) > 0.01) | (np.abs(mlow-np.median(mlow)) > 0.015)
    bad = (np.abs(mmid-np.median(mmid)) > 0.1)
    xx = np.arange(len(mmid))
    s = np.argsort(xglint[~bad,1])
    ds9.view((yarr.T / mmid)[:,~bad][:,s].T)
    
    s = np.argsort(yglint[~bad])
    
    #######################
    ####### Figure for ISR
    ###
    label = 'G141, zodi'
    #label = 'G102, zodi'
    label = 'G141, He excess'
    bad = yglint > 0.98
    s = np.argsort(xglint[~bad,1])
    
    label = 'G102, He excess'
    
    #plt.gray()
    fig = unicorn.plotting.plot_init(aspect=1, square=True, hspace=0, left=0.1, bottom=0.1, use_tex=True, NO_GUI=True, xs=5, fontsize=11)
    
    dx, left = 0.57, 0.115
    ax = fig.add_axes((left, 0.1, dx, 0.88))
    ax.imshow(nd.median_filter((yarr.T / mmid)[:,~bad][:,s].T, size=4), aspect='auto', vmin=0.85, vmax=1.1)
    ax.set_ylabel(r'$N_\mathrm{read}$')
    ax.set_xlabel(r'Image column, $x$')
    
    ddx = 0.03
    ax2 = fig.add_axes((left+dx+ddx, 0.1, (0.98-dx-left-ddx), 0.88))
    n = np.arange((~bad).sum())
    ax2.plot(xglint[~bad,1][np.argsort(xglint[~bad,1])], n, color='black')
    ax2.set_ylim(0, (~bad).sum())
    ax2.set_xlim(0.3, 3.05)
    ax2.set_xlabel(r'Background (e$^{-}$/s)')
    ax2.set_yticklabels([])
    ax.text(0.05, 0.95, label, ha='left', va='top', backgroundcolor='white', bbox=dict(facecolor='white', edgecolor='None', alpha=0.5), transform=ax.transAxes, size=12)
    
    ax2.set_xlim(0.2, 5)
    
    unicorn.plotting.savefig(fig, '/tmp/column_average_%s.pdf' %(label.replace(', ', '_').replace(' ','_')))
    
    
    #         
    # for file in files[::(len(files)/200)]:
    #     print file
    #     medx = pyfits.open(file)
    #     med = medx[0]
    #     #plt.scatter(xarr, med.data, marker='.', c=med.header['BGOBS']*(xarr*0+1), vmin=0.5, vmax=2, alpha=0.1)
    #     plt.plot(xarr, med.data/np.median(med.data[-150:]), alpha=0.01, color='black') #, color='%.3f' %(cm(med.header['BGOBS'])), alpha=0.1)
    #     medx.close()
        
    ### Make medians on subimages 
    root, glint_limit = 'low', 1.1
    root, glint_limit = 'high', 1.02 #0.98

    files=glob.glob('i*%s_??.fits' %(root)) ### All
    #files=np.array(files)[yglint < 0.98]
    files=np.array(files)[(yglint > 0.98) & (yglint < 1.3)]
    #files=glob.glob('i?[a-z]*%s_??.fits' %(root)) ### No GOODS-N
    
    files=glob.glob('ib*%s_??.fits' %(root)) ### GOODS-N high
    #files=glob.glob('ibh*%s_??.fits' %(root)) ### 3D-HST high
    files=glob.glob('ibhm[3-9]*%s_??.fits' %(root)) ### COSMOS low
    files=glob.glob('ibhj*%s_??.fits' %(root)) ### GOODS-S/AEGIS
    #
    #files=glob.glob('ic*%s_??.fits' %(root)) ### G102
    #files=glob.glob('ic*%s_??.fits' %(root)) ### G102
    
    files = glob.glob('i*%s_??.fits' %(root)) ### G102, low
    files = glob.glob('ibh*%s_??.fits' %(root)) ### G102, low
    #files = glob.glob('ic14*%s_??.fits' %(root)) ### G102, low
    
    ### Check the images
    base=''
    for file in files:
        base_i = file.split('_')[0]
        if base_i == base:
            continue
        #
        base = base_i
        im = pyfits.open(file)
        ok = im[0].data > -90
        med = np.median(im[0].data[ok])
        ds9.view(im[0].data-med)
        ds9.scale(-1,2)
        if os.path.exists('%s_addmask.reg' %(base)):
            ds9.set('regions file %s_addmask.reg' %(base))
        #    
        xx = raw_input('%s, mask_ok? ' %(file))
        ds9.set('regions save %s_addmask.reg' %(base))
        if xx == 'q':
            break
        #
        if xx != '':
            os.system('echo "mv %s*%s* tmp/" >> %s_bad.list' %(base, root, root))
    
    #### Apply new masks
    masks = glob.glob('*addmask.reg')
    for mask in masks:
        lines = open(mask).readlines()
        if len(lines) < 4:
            continue
        #
        base = mask.split('_ad')[0]
        images = glob.glob('%s_%s_??.fits' %(base, root))
        if len(images) == 0:
            continue
        #
        im = pyfits.open(images[0])
        add_mask = pyregion.open('%s_addmask.reg' %(base)).as_imagecoord(header=im[0].header).get_mask(hdu=im[0])
        if add_mask.sum() == 0:
            continue
        #
        for image in images:
            print 'Fix image %s' %(image)
            im = pyfits.open(image, mode='update')
            im[0].data[add_mask > 0] = -100
            im.flush()
    
    ##### Combine into master image       
    master = np.zeros((1014,1014))
    master_wht = np.zeros((1014,1014))
    master_std = np.zeros((1014,1014))
    Npix = 169 # 100/200
    big = np.zeros((len(files), Npix, Npix))
    for i in range(1014/Npix):
        for j in range(1014/Npix):
            sly, slx = slice(i*Npix,(i+1)*Npix), slice(j*Npix,(j+1)*Npix)
            for k in range(len(files)):
                im = pyfits.open(files[k])
                subim = im[0].data[sly, slx]
                ok = (subim > -90) & (subim < 5)
                nmad = threedhst.utils.nmad(subim[ok])
                med = np.median(subim[ok])
                bad = np.abs(subim-med) > 4*nmad
                subim[bad] = -100              
                big[k,:,:] = subim
                med_file = pyfits.open(files[k].replace('.fits','_med.fits'))
                med = med_file[0].data
                #whts[k] = med_file[0].header['BGSCL']
                glint = np.median(med[200:400])/np.median(med[750:900])
                mright = np.median(med[500:900])
                subim[ok] /= mright
                big[k,:,:] = subim
                if glint > glint_limit:
                    big[k,:,:] = -100
                    print 'glint'
                #
                print '(%d,%d): %04d, %s' %(i,j,k, files[k])
            #
            mask = (big < -90) | (big > 3)
            limits = np.percentile(big[~mask].flatten(), [2.5, 97.5])
            mask &= (big < limits[0]) | (big > limits[1])
            mx = ma.masked_array(big, mask=mask, fill_value=-100)
            #
            med = ma.mean(mx, axis=0)
            master[sly, slx] = med*1
            pyfits.writeto('master_%s.fits' %(root), data=master, clobber=True)
            #master_wht[sly, slx] = ((big > -90).sum(axis=0))
            master_wht[sly, slx] = len(files)-mx.mask.sum(axis=0)
            pyfits.writeto('master_%s_wht.fits' %(root), data=master_wht, clobber=True)
            master_std[sly, slx] = ma.std(mx, axis=0)
            pyfits.writeto('master_%s_std.fits' %(root), data=master_std, clobber=True)
    
    ### last strips
    strip = 1014-1014/Npix*Npix
    big = np.zeros((len(files), strip, 1014))
    sly, slx = slice(1014-strip, 1014), slice(0,1014)
    ### copy from "for k..." above
    
    big = np.zeros((len(files), 1014, strip))
    slx, sly = slice(1014-strip, 1014), slice(0,1014)
    ### copy from "for k..." above
    
    ##
    pyfits.writeto('zodi_G102.fits', data=master, clobber=True)
    pyfits.writeto('zodi_G102_wht.fits', data=master_wht, clobber=True)
    
    pyfits.writeto('excess_G102.fits', data=master, clobber=True)
    pyfits.writeto('excess_G102_wht.fits', data=master_wht, clobber=True)

    pyfits.writeto('zodi_G141.fits', data=master, clobber=True)
    pyfits.writeto('zodi_G141_wht.fits', data=master_wht, clobber=True)

    pyfits.writeto('excess_G141.fits', data=master, clobber=True)
    pyfits.writeto('excess_G141_wht.fits', data=master_wht, clobber=True)
    pyfits.writeto('excess_G141_std.fits', data=master_std, clobber=True)

    pyfits.writeto('excess_lo_G141.fits', data=master, clobber=True)
    pyfits.writeto('excess_lo_G141_wht.fits', data=master_wht, clobber=True)
    pyfits.writeto('excess_lo_G141_std.fits', data=master_std, clobber=True)

    pyfits.writeto('excess_sl_G141.fits', data=master, clobber=True)
    pyfits.writeto('excess_sl_G141_wht.fits', data=master_wht, clobber=True)
    pyfits.writeto('excess_sl_G141_std.fits', data=master_std, clobber=True)
    
    pyfits.writeto('excess_3dhst_G141.fits', data=master, clobber=True)
    pyfits.writeto('excess_3dhst_G141_wht.fits', data=master_wht, clobber=True)

    pyfits.writeto('excess_goodsn_G141.fits', data=master, clobber=True)
    pyfits.writeto('excess_goodsn_G141_wht.fits', data=master_wht, clobber=True)
    
    #### scattered light
    ex = pyfits.open('excess_lo_G141.fits')[0].data
    sl = pyfits.open('excess_sl_G141.fits')[0].data
    diff = m-ex
    sy = np.median(diff[:,50:200], axis=1)
    sx = np.median(((m-ex).T/sy).T, axis=0)
    xarr = np.arange(1014)
    from scipy import polyfit, polyval
    cx = polyfit(xarr-507, sx, 5)
    cy = polyfit(xarr-507, sy, 5)
    model_sl = ((np.ones((1014, 1014)) * polyval(cx, xarr-507)).T*(polyval(cy, xarr-507))).T
    pyfits.writeto('G141_scattered_light.fits', data=model_sl, clobber=True)
    
    #### Clean up final version
    flat_140 = pyfits.open('%s/uc721143i_pfl.fits' %(os.getenv('iref')))[1].data
    flat_g141 = pyfits.open('%s/u4m1335mi_pfl.fits' %(os.getenv('iref')))[1].data
    flat = flat_140/flat_g141

    flat_105 = pyfits.open('%s/uc72113oi_pfl.fits' %(os.getenv('iref')))[1].data
    flat_125 = pyfits.open('%s/uc72113qi_pfl.fits' %(os.getenv('iref')))[1].data
    flat_160 = pyfits.open('%s/uc721145i_pfl.fits' %(os.getenv('iref')))[1].data
    f160 = pyfits.open('%s/uc721145i_pfl.fits' %(os.getenv('iref')))
    
    better_flat = (0.3*flat_125+0.7*flat_160)
    f160[1].data = better_flat
    f160[0].header.add_comment('Zodi optimized flat')
    f160[0].header.add_comment('G. Brammer (%s)' %(time.ctime()))
    f160[0].header.add_comment('0.3*F125W + 0.7*F160W')
    f160[0].header.add_comment('F125W: uc72113qi_pfl')
    f160[0].header.add_comment('F160W: uc721145i_pfl')
    f160.writeto('zodi_G141_flat.fits', clobber=True)

    better_flat_he = (0.7*flat_105+0.3*flat_160)
    f160[1].data = better_flat_he
    f160.writeto('excess_G141_flat.fits', clobber=True)

    import time
    f105 = pyfits.open('%s/uc72113oi_pfl.fits' %(os.getenv('iref')))
    better_flat_he = (0.7*flat_105+0.3*flat_125)
    f105[1].data = better_flat_he
    f105[0].header.add_comment('He10830 optimized flat')
    f105[0].header.add_comment('G. Brammer (%s)' %(time.ctime()))
    f105[0].header.add_comment('0.7*F105W + 0.3*F125W')
    f105[0].header.add_comment('F105W: uc72113oi_pfl')
    f105[0].header.add_comment('F125W: uc72113qi_pfl')
    
    f105.writeto('excess_G102_flat.fits', clobber=True)
    
    ref = pyfits.open('%s/CONF/WFC3.IR.G141.sky.V1.0.fits' %(os.getenv('THREEDHST')))   
    
    new = pyfits.open('zodi_G141.fits')
    new_wht = pyfits.open('zodi_G141_wht.fits')

    new = pyfits.open('excess_G141.fits')
    new_wht = pyfits.open('excess_G141_wht.fits')

    new = pyfits.open('excess_lo_G141.fits')
    new_wht = pyfits.open('excess_lo_G141_wht.fits')

    new = pyfits.open('zodi_G102.fits')
    new_wht = pyfits.open('zodi_G102_wht.fits')

    new = pyfits.open('excess_G102.fits')
    new_wht = pyfits.open('excess_G102_wht.fits')

    new[0].data *= (flat_140/better_flat)[5:-5, 5:-5]
    f160[1].data = better_flat
    f160.writeto('zodi_G141_flat.fits', clobber=True)

    new = pyfits.open('excess_G141.fits')
    new_wht = pyfits.open('excess_G141_wht.fits')

    new = pyfits.open('excess_goodsn_G141.fits')
    new_wht = pyfits.open('excess_goodsn_G141_wht.fits')

    new = pyfits.open('excess2_G141.fits')
    new_wht = pyfits.open('excess2_G141_wht.fits')
    
    new[0].data *= (better_flat/better_flat_he)[5:-5, 5:-5]
    
    for i in range(5):
        med = nd.median_filter(new[0].data, size=8)
        #fill = new_wht[0].data == 0
        kernel = np.ones((7,7))
        sigma = np.sqrt(nd.convolve((new[0].data-med)**2, kernel)/kernel.sum())
        #
        bad = (new_wht[0].data == 0) | (np.abs(new[0].data/med-1) > 4*sigma)
        new[0].data[bad] = med[bad]
            
    flt = pyfits.open('/Users/brammer/WFC3/Backgrounds/MaskedMaster/G141/ibhj01ioq_flt.fits.gz')
    deathstar = (flt['DQ'].data & 4) > 0
    new[0].data[deathstar] = 1
    new[0].data[new[0].data == 0] = med[new[0].data == 0]
    new[0].data[new[0].data == 0] = 1
    new_wht[0].data[deathstar] = 0
    
    dq = (bad | deathstar)*1
    
    # root='zodi'
    # root='excess'
    # root='excess2'
    root = new.filename().split('.fits')[0]
        
    new.writeto('%s_clean.fits' %(root), clobber=True)
    new_wht.writeto('%s_clean_wht.fits' %(root), clobber=True)
    pyfits.writeto('%s_clean_dq.fits' %(root), data=dq, clobber=True)
    
    GRISM = 'G102'
    GRISM = 'G141'
    
    ims = [pyfits.open('zodi_%s_clean.fits' %(GRISM))[0].data, pyfits.open('excess_%s_clean.fits' %(GRISM))[0].data]
    for i in range(2):
        ims[i][ims[i] == 0] = 1
    
    for f in np.arange(0,1.01, 0.25):
        sum = ims[1]*f + ims[0]*(1-f)
        pyfits.writeto('%s_zodi_%4.2f_excess_%4.2f.fits' %(GRISM, 1-f, f), clobber=True, data=sum)
        
    ##########
    
    
    #### Image statistics
    xx = pyfits.open('/Users/brammer/3DHST/Spectra/Work/CONF/WFC3.IR.G141.sky.V1.0.fits')   
    yy = pyfits.open('/Users/brammer/3DHST/Spectra/Work/CONF/sky.G141.set001.fits')   
    
    slx, sly = slice(670, 770), slice(103, 203)
    ref_im = (xx[0].data*(1/better_flat*flat_g141)[5:-5,5:-5])[sly, slx].flatten()
    my_im = (yy[0].data)[sly, slx].flatten()
    new_im = (low[0].data*(flat_140/better_flat)[5:-5,5:-5])[sly, slx].flatten()
    #ref_im = (xx[0].data*(flat_g141)[5:-5,5:-5])[sly, slx].flatten()
    #new_im = (low[0].data*(flat_140)[5:-5,5:-5])[sly, slx].flatten()
    
    #flat_im = better_flat[5:-5,5:-5][sly, slx].flatten()
    
    plt.hist(ref_im/np.median(ref_im), bins=100, range=(0.96,1.04), alpha=0.8)
    plt.hist(my_im/np.median(my_im), bins=100, range=(0.96,1.04), alpha=0.8)
    plt.hist(new_im/np.median(new_im), bins=100, range=(0.96,1.04), alpha=0.8)
    #plt.hist(flat_im/np.median(flat_im), bins=100, range=(0.9,1.1), alpha=0.8)
    
    print np.diff(np.percentile(ref_im, [16,84]))/2
    print np.diff(np.percentile(my_im, [16,84]))/2
    print np.diff(np.percentile(new_im, [16,84]))/2
    #print np.diff(np.percentile(flat_im, [16,84]))/2
    
    #### Check correction
    zodi = pyfits.open('zodi_G141_clean.fits')[0].data
    zodi_dq = pyfits.open('zodi_G141_clean_dq.fits')[0].data
    zodi_flat = pyfits.open('zodi_G141_flat.fits')[1].data[5:-5,5:-5]
    excess = pyfits.open('excess_G141_clean.fits')[0].data
    excess_dq = pyfits.open('excess_G141_clean_dq.fits')[0].data
    excess_flat = pyfits.open('excess_G141_flat.fits')[1].data[5:-5,5:-5]
    excess2 = pyfits.open('excess2_G141_clean.fits')[0].data
    excess2_dq = pyfits.open('excess2_G141_clean_dq.fits')[0].data
    excess2_flat = pyfits.open('excess_G141_flat.fits')[1].data[5:-5,5:-5]
    
    #ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj02l0q_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj44mrq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj47cfq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj47cmq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj49aiq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj38a6q_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj38adq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj31giq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj31grq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ib3720f0q_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ib3720f4q_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ib3715s7q_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ib3715sbq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ib3708i9q_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ib3706b6q_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj70wlq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj42o5q_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj42ocq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhm19pqq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhm19pxq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhm30zsq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhm42umq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhm46iaq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhm54doq_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj20x0q_ima.fits')
    ima = pyfits.open('/Volumes/3DHST_Gabe/MasterBackground/IMA/ibhj20x7q_ima.fits')
    
    bg = catIO.Readfile('../RAW/%sj_G141_orbit.dat' %(os.path.basename(ima.filename().split('q_ima')[0])))
    
    mask_file = '/Volumes/3DHST_Gabe/MasterBackground/IMA/%sq_flt_mask.fits' %(os.path.basename(ima.filename().split('q_ima')[0]))
    mask = (pyfits.open(mask_file)[0].data > 0) | (zodi == 0) | (excess == 0) | (zodi_dq > 0) | (excess_dq > 0) | (excess2 == 0)
    
    sci = np.cast[np.double](ima['sci',1].data[5:-5,5:-5][~mask].flatten())
    var = np.cast[np.double](ima['err',1].data[5:-5,5:-5][~mask].flatten())**2
    
    templates = np.array([(zodi*zodi_flat).flatten(), (excess*excess_flat).flatten()])
    flats = np.array([zodi_flat.flatten(), excess_flat.flatten()])
    init = [bg.zodi[0], np.maximum(np.mean(bg.bg-bg.zodi), 0.02)]
    
    templates = np.array([(zodi*zodi_flat).flatten(), (excess*excess_flat).flatten(), (excess2*excess_flat).flatten()])
    flats = np.array([zodi_flat.flatten(), excess_flat.flatten(), excess_flat.flatten()])
    init = [bg.zodi[0], np.maximum(np.mean(bg.bg-bg.zodi), 0.02)/2., np.maximum(np.mean(bg.bg-bg.zodi), 0.02)/2.]
    
    # templates = np.array([(zodi*zodi_flat).flatten(), (excess2*excess_flat).flatten()])
    # flats = np.array([zodi_flat.flatten(), excess_flat.flatten()])
    # init = [bg.zodi[0], np.maximum(np.mean(bg.bg-bg.zodi), 0.02)]
    
    amatrix = unicorn.utils_c.prepare_nmf_amatrix(var, templates[:,(~mask).flatten()])
    coeffs = unicorn.utils_c.run_nmf(sci, var, templates[:,(~mask).flatten()], amatrix, toler=1.e-6, MAXITER=100000, verbose=1, init_coeffs=init)
        
    model_flat = np.dot(coeffs.reshape((1,-1)), flats).reshape((1014,1014))/coeffs.sum()
    
    model = np.dot(coeffs.reshape((1,-1)), templates).reshape((1014,1014))/model_flat
    
    med = np.median(model)
    
    ds9.frame(1)
    ds9.view((ima['sci',1].data*flat_g141)[5:-5,5:-5]/model_flat*(~mask))
    ds9.scale(med-0.5, med+0.5)

    ds9.frame(2)
    ds9.view(model)
    ds9.scale(med-0.5, med+0.5)

    ds9.frame(3)
    ds9.view((ima['sci',1].data*flat_g141)[5:-5,5:-5]/model_flat*(~mask) - model*(~mask))
    ds9.scale(-0.5, 0.5)
    
    flt = pyfits.open('/Users/brammer/3DHST/Spectra/Work/3DHST_VariableBackgrounds/GOODS-S/%s' %(os.path.basename(ima.filename()).replace('ima','flt')))[1].data
    ds9.frame(4)
    ds9.view(flt*(~mask))
    ds9.scale(-0.5, 0.5)
    
def visit_minimum(base='ic1409y2q', PATH='./IMA/'):
    import os
    from threedhst import catIO
    
    visit = os.path.basename(base)[:6]
    if not os.path.exists('%s_visit.dat' %(visit)):
        os.system('cat %s/%s*orbit.dat > %s_visit.dat' %(PATH, visit, visit))
    
    vs = catIO.Readfile('%s_visit.dat' %(visit))
    ok = (vs.read != 1) & (vs.shadow == 1) & (vs.limbang > 30)
    
    if ok.sum() == 0:
        zodi = vs.zodi.min()
    else:
        zodi = vs.bg[ok].min()
    
    return zodi, vs.bg[ok]
    
def remove_multi_sky(flt='icat21dgq_flt.fits', list=['zodi_G102_clean.fits', 'excess_G102_clean.fits'],  path_to_sky = '../CONF/', out_path='./', verbose=False, plot=False, flat_correct=True, sky_subtract=True, second_pass=True, overall=True, combine_skies=False, grow_mask=10):
    import scipy.ndimage as nd
    import threedhst.grism_sky as bg
    #import scipy.signal as sign
    
    # flt = '../../GOODS-N/RAW/ib3708ilq_flt.fits.gz'
    im = pyfits.open(flt)
    bg.set_grism_flat(grism=im[0].header['FILTER'])
    
    segfile = os.path.basename(flt.replace('.fits','.seg.fits')).replace('.gz','')
    if os.path.exists(segfile):
        seg = pyfits.open(segfile)[0].data
        use_biweight=False
    else:
        seg = np.zeros(im[1].data.shape)
        use_biweight=True
    
    mask = nd.maximum_filter((seg > 0)*1, size=grow_mask)
    
    templates = []
    for i in range(len(list)):
        templates.append(pyfits.open(path_to_sky + list[i])[0].data.flatten())
        
    ok = (mask == 0) & ((im['DQ'].data & (4+16+32)) == 0) & (im['ERR'].data > 0) & (np.sum(templates, axis=0) > 0).reshape((1014, 1014))
    
    sn_med = np.median((im['SCI'].data/im['ERR'].data)[ok])
    ok &= im['SCI'].data/im['ERR'].data < 2*sn_med
    
    sci = np.cast[np.double]((im['SCI'].data*bg.flat)[ok]).flatten()
    var = np.cast[np.double](im['ERR'].data[ok]).flatten()**2
    templates = np.cast[np.double](templates)#[:,ok.flatten()]
    
    NTEMP = len(list)
    
    init = np.ones(NTEMP)*np.median(sci)/NTEMP
    
    amatrix = unicorn.utils_c.prepare_nmf_amatrix(var, templates[:, ok.flatten()])    
    coeffs = unicorn.utils_c.run_nmf(sci, var, templates[:, ok.flatten()], amatrix, toler=1.e-7, MAXITER=100000, verbose=1, init_coeffs=init)
    
    model = np.dot(coeffs, templates).reshape((1014, 1014))
    
    pyfits.writeto('xx.fits', data=model, clobber=True)
    pyfits.writeto('yy.fits', data=(im['SCI'].data*bg.flat), clobber=True)
    pyfits.writeto('zz.fits', data=(im['SCI'].data*bg.flat)-model, clobber=True)
    
    im['SCI'].data *= bg.flat
    im.writeto('yy_flt.fits', clobber=True)
    xin, yin = bg.profile('yy_flt.fits', extension=1, flatcorr=False, biweight=use_biweight)

    im['SCI'].data -= model
    im.writeto('zz_flt.fits', clobber=True)
    x2, y2 = bg.profile('zz_flt.fits', extension=1, flatcorr=False, biweight=use_biweight)
    
    im['SCI'].data = model*1
    im.writeto('mm_flt.fits', clobber=True)
    xm, ym = bg.profile('mm_flt.fits', extension=1, flatcorr=False, biweight=use_biweight)
    
    
    
def mask_visit(root='AEGIS-10', threshold=0.005, grow_mask=8):
    """
    Extract masks from 3D-HST model images and blot them 
    back to the FLT frames of the individual exposures.
    
    """
    import scipy.ndimage as nd
    import pyfits
    import threedhst
    import unicorn
    import numpy as np
    
    # root='AEGIS-10'; threshold=0.005; grow_mask=8
    
    asn = threedhst.utils.ASNFile('%s-G141_asn.fits' %(root.lower().replace('goods-','goods')))
    model = pyfits.open('%s_inter_model.fits' %(root))
    
    xoff, yoff = unicorn.reduce.get_interlace_offsets('%s-G141_asn.fits' %(root), growx=model[0].header['GROWX'], growy=model[0].header['GROWY'])
    
    mask = np.cast[np.uint8]((nd.maximum_filter((model[0].data > threshold)*1, size=grow_mask) > 0)*1)
    #mask = model[0].data
    
    for i in range(len(asn.exposures)): 
        #flt = pyfits.open('%s_flt.fits.gz' %(asn.exposures[i]))
        x0 = model[0].header['NGROW']*model[0].header['GROWX'] + model[0].header['PAD']/2 + xoff[i]
        y0 = model[0].header['NGROW']*model[0].header['GROWY'] + model[0].header['PAD']/2 + yoff[i]
        sub_mask = mask[y0:y0+1014*model[0].header['GROWY']:model[0].header['GROWY'], x0:x0+1014*model[0].header['GROWX']:model[0].header['GROWX']]
        pyfits.writeto('%s_flt_mask.fits' %(asn.exposures[i]), data=sub_mask, clobber=True)
        print asn.exposures[i]
        
def simple_mask():
    """
    Try making a grism mask with a simple convolution kernel
    """
    import os
    import pyfits
    import unicorn
    import stsci.image 
    
    os.chdir("/Users/brammer/3DHST/Spectra/Work/Rudnick/PREP_FLT")
    
    direct = pyfits.open('ic1b04o7q_flt.fits')
    grism = pyfits.open('ic1b04o9q_flt.fits')
    grism_seg = pyfits.open('ic1b04o9q_flt.seg.fits')
    
    #
    unicorn.reduce.set_grism_config(grism='G102', force=True)
    orders, xi = unicorn.reduce.grism_model(507, 507, lam_spec=None, flux_spec=None, BEAMS=['A', 'B'], grow_factor=1, growx=1, growy=1, pad=0, ngrow=0, grism='G102')
    
    yord, xord = np.indices(orders.shape)
    beams = np.dot(np.ones((orders.shape[0],1), dtype=np.int), xi[5].reshape((1,-1)))
    
    cast = np.int
    xord, yord, ford, word, sord, bord = np.array(xord, dtype=cast), np.array(yord, dtype=cast), np.array(orders, dtype=np.float64), xi[2], xi[3], beams
    
    data = direct[1].data*1
    
    mask = (data > 10) | (data < -10) | ((direct['DQ'].data & (4+16+32)) > 0)
    
    data[mask] = 0
    bg = np.median(data[~mask])
    
    SN = (direct['SCI'].data - bg) / direct['ERR'].data
    mask |= (SN < 0.2) 
    
    data[mask] = 0
    bg = np.median(data[~mask])
    data[mask] = bg*1
    
    sh = data.shape
    dpad = 300
    pad = np.zeros((sh[0], sh[1]+2*dpad))
    slx = slice(dpad, -dpad)
    pad[:,slx] += data-bg
    
    model = stsci.image.convolve.convolve2d(pad, ford/254, fft=True)
    dx = ford.shape[1]/2 + xi[0]
    cut_model = model[:, dpad-dx:dpad-dx+1014]
    
    grism_mask = (cut_model > 0.02)*1.
    NGROW = 12
    grism_mask_grow = nd.maximum_filter(grism_mask, size=NGROW)
    pyfits.writeto(grism.filename().replace('flt','flt.mask'), data=np.cast[np.uint8](grism_mask_grow), header=grism[1].header, clobber=True)
    
    # non_zero = orders > 0
    # non_zero = orders > -100
    # 
    # cast = np.int
    # xord, yord, ford, word, sord, bord = np.array(xord[non_zero], dtype=cast), np.array(yord[non_zero], dtype=cast), np.array(orders[non_zero], dtype=np.float64), xi[2][non_zero], xi[3][non_zero], beams[non_zero]
    
######## Figures for the ISR

def track_figure(visit='ibhj07', PATH='/Volumes/3DHST_Gabe/MasterBackground'):
    """
    Make a figure showing which reads satisfy the criteria for the master 
    images.
    """
    
    bg = catIO.Readfile('%s/%s_visit.dat' %(PATH, visit))
    times = bg.seconds*0.
    for i in range(bg.N):
        ima = pyfits.open('%s/IMA/%sq_ima.fits' %(PATH, bg.name[i].split('j_')[0]))
        times[i] = bg.seconds[i]/86400 + ima[0].header['EXPSTART']
        print i, bg.N
    
    fig = unicorn.plotting.plot_init(aspect=1/3., xs=7, fontsize=11, left=0.1, bottom=0.084, use_tex=True, NO_GUI=True, right=0.03)
    ax = fig.add_subplot(111)
    
    t0 = (times-times[0])*86400/60.
    ix = np.arange(bg.N)[bg.read == 1]
    ix = np.append(ix, bg.N+1)
    for i in range(len(ix)-1):
        sl = slice(ix[i], ix[i+1])
        ax.plot(t0[sl], bg.bg[sl], color='black', alpha=0.5, linewidth=5, zorder=-1)
    
    ax.scatter(t0, bg.bg, color='white', s=60, zorder=1)
    ax.scatter(t0, bg.bg, color='black', s=40, alpha=0.8, label='SPARS100 reads')
    zodi = bg.shadow == 1
    ax.scatter(t0[zodi], bg.bg[zodi], color='green', s=40, alpha=0.8, label='Shadow (zodi-dominated)')
    
    ax.plot(t0, bg.zodi, color='blue', linestyle='--', linewidth=2, label='ETC zodi model')
    
    zodi_level, zodi_all = mywfc3.masterbg.visit_minimum(base=visit, PATH='./IMA/')
    level, ok_reads = 'high', np.arange(bg.N)[(bg.bg-zodi_level > 0.7) & (bg.shadow == 0)] + 1
    ax.scatter(t0[ok_reads-1], bg.bg[ok_reads-1], color='red', s=40, alpha=0.8, label=r'He 10830\,\AA\ $> 0.7$ e$^-$/s')
    
    ax.set_xlim(-5, 140)
    ax.set_ylim(0,3.1)
    ax.set_xlabel(r'$\Delta t$ (minutes)')
    ax.set_ylabel(r'Background (e$^{-}$/s)')
    ax.legend(loc=(0.4,0.4), scatterpoints=1, title='Visit: %s (G141)' %(visit), prop=dict(size=9))
    
    unicorn.plotting.savefig(fig, '/tmp/masterbg_track.pdf')
    
    
def example_exposure(root='ibhj07ynq'):
    """
    Show an example
    """
    ima = pyfits.open('IMA/%s_ima.fits' %(root))['sci',1].data[5:-5,5:-5]
    
    flat_f140 = pyfits.open('%s/uc721143i_pfl.fits' %(os.getenv('iref')))[1].data[5:-5,5:-5]
    flat_g141 = pyfits.open('%s/u4m1335mi_pfl.fits' %(os.getenv('iref')))[1].data[5:-5,5:-5]
    
    mask = pyfits.open('IMA/%s_flt_mask.fits' %(root))[0].data
    
    med = np.median((ima/flat_f140)[mask == 0])
    
    fig = unicorn.plotting.plot_init(aspect=1/3., xs=7, fontsize=11, left=0.01, bottom=0.01, use_tex=True, NO_GUI=True, right=0.01, top=0.01, wspace=0.02, hspace=0.02)
    
    v = [0.85, 1.1]
    
    ax = fig.add_subplot(131)
    ax.imshow(ima/med, vmin=v[0], vmax=v[1], aspect='auto', interpolation='nearest')
    ax.text(0.02, 0.98, 'a) ' + root + ' (G141)', ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', pad=6), fontsize=10,  transform=ax.transAxes)
    
    ax = fig.add_subplot(132)
    ax.imshow(ima/med/flat_f140, vmin=v[0], vmax=v[1], aspect='auto', interpolation='nearest')
    ax.text(0.02, 0.98, 'b) Divide by F140W flat', ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', pad=6), fontsize=10,  transform=ax.transAxes)

    ax = fig.add_subplot(133)
    ax.imshow(ima/med/flat_f140*(mask == 0), vmin=v[0], vmax=v[1], aspect='auto', interpolation='nearest')
    ax.text(0.02, 0.98, 'c) Object mask', ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', pad=6), fontsize=10,  transform=ax.transAxes)
    
    for ax in fig.axes:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    
    unicorn.plotting.savefig(fig, '/tmp/image_example.pdf', dpi=200)
    
def show_master_images():
    
    v = [0.75, 1.08]
    v = [0.8,1.1]
    
    GRISM = 'G141'    
    zodi = pyfits.open('zodi_G141_clean.fits')[0].data
    zodi_wht = pyfits.open('zodi_G141_clean_wht.fits')[0].data

    excess = pyfits.open('excess_lo_G141_clean.fits')[0].data
    excess_wht = pyfits.open('excess_lo_G141_clean_wht.fits')[0].data
    

    GRISM = 'G102'
    zodi = pyfits.open('../MasterBackgroundG102/zodi_G102_clean.fits')[0].data
    zodi_wht = pyfits.open('../MasterBackgroundG102/zodi_G102_clean_wht.fits')[0].data
    
    excess = pyfits.open('../MasterBackgroundG102/excess_G102_clean.fits')[0].data
    excess_wht = pyfits.open('../MasterBackgroundG102/excess_G102_clean_wht.fits')[0].data
    # 
    flt = pyfits.open('/Users/brammer/WFC3/Backgrounds/MaskedMaster/G141/ibhj01ioq_flt.fits.gz')
    deathstar = (flt['DQ'].data & 4) == 0
    
    hist_pad = 0.3
    
    #v = [0.85, 1.1]
    
    fig = unicorn.plotting.plot_init(aspect=1/2.*(1+hist_pad), xs=7, fontsize=11, left=0.01, bottom=0.01, use_tex=True, NO_GUI=True, right=0.01, top=0.01, wspace=0.02, hspace=0.02)
    
    left, ddx, bottom, right, top = 0.01, 0.01, 0.085, 0.9, 0.98
    dx = (right-left-ddx)/2.    
    pp = [2,98]
    
    ax = fig.add_axes((left, hist_pad, dx, top-hist_pad))
    ax.imshow(zodi*deathstar, vmin=v[0], vmax=v[1], aspect='auto', interpolation='nearest')
    ax.set_xticklabels([]); ax.set_yticklabels([])

    ax = fig.add_axes((left, bottom, dx, hist_pad-bottom-0.03))
    ok = (zodi_wht != 1) & (zodi_wht > 0)
    limits = np.percentile(zodi_wht[ok], pp)
    ax.hist(zodi_wht[ok].flatten(), bins=np.minimum(50, limits[1]-limits[0]), alpha=0.5, color='black', range=limits, histtype='stepfilled')
    ax.set_xlim(limits)
    ax.set_yticklabels([])
    ax.set_xlabel('Reads / pixel')
    ax.text(0.03, 0.92, '%s - zodi' %(GRISM), va='top', ha='left', transform=ax.transAxes)
    
    ax = fig.add_axes((left+dx+ddx, hist_pad, dx, top-hist_pad))
    imsh = ax.imshow(excess*deathstar, vmin=v[0], vmax=v[1], aspect='auto', interpolation='nearest')
    ax.set_xticklabels([]); ax.set_yticklabels([])
    cax = fig.add_axes((left+dx+ddx+dx+ddx, hist_pad, 0.02, top-hist_pad))
    cb = fig.colorbar(imsh, cax=cax)
    
    ax = fig.add_axes((left+dx+ddx, bottom, dx, hist_pad-bottom-0.03))
    ok = (zodi_wht != 1) & (zodi_wht > 0)
    limits = np.percentile(excess_wht[ok], pp)
    ax.hist(excess_wht[ok].flatten(), bins=np.minimum(50, limits[1]-limits[0]), alpha=0.5, color='black', range=limits, histtype='stepfilled')
    ax.set_xlim(limits)
    ax.set_yticklabels([])
    ax.set_xlabel('Reads / pixel')
    ax.text(0.03, 0.92, '%s - He' %(GRISM), va='top', ha='left', transform=ax.transAxes)
    
    unicorn.plotting.savefig(fig, '/tmp/master_image_%s.pdf' %(GRISM), dpi=300)
    
def background_components():
    """
    Show how different orders contribute to the grism background
    """
    import pyfits
    import numpy as np
    
    import pysynphot as S
    import unicorn
    import matplotlib.pyplot as plt
    
    grism = 'G141'
    ex = pyfits.open('excess_lo_G141_clean.fits')[0].data
    #ex = pyfits.open('zodi_G141_clean.fits')[0].data
    med_ex = np.median(ex[480:520,:], axis=0)
    
    # grism = 'G102'
    # ex2 = pyfits.open('../MasterBackgroundG102/excess_G102_clean.fits')[0].data
    # ex2 = pyfits.open('../MasterBackgroundG102/zodi_G102_clean.fits')[0].data
    # med_ex = np.median(ex2[480:520,:], axis=0)
    
    heline = S.GaussianSource(1.e-16, 1.083e4, 10)
    #heline = S.GaussianSource(1.e-16, 1.283e4, 10)
    continuum = S.FlatSpectrum(29.9, fluxunits='STMag')
    #continuum = S.FlatSpectrum(23, fluxunits='ABMag')
    #x = np.arange(5000,1.8e4,2)
    #beta, n = -1.5, 1.e-18
    #continuum = S.ArraySpectrum(x, (x/1.e4)**beta*n, fluxunits='flam')
    obs = heline+continuum
    
    #obs = S.FileSpectrum('$PYSYN_CDBS/etc/background/zodiacal_model_001.fits')
    
    unicorn.reduce.set_grism_config(grism=grism)
    orders, xi = unicorn.reduce.grism_model(507, 507, lam_spec=np.cast[float](obs.wave), flux_spec=np.cast[float](obs.flux)/1.e-17, BEAMS=['A', 'B', 'C', 'D', 'E'], grow_factor=1, growx=1, growy=1, pad=0, ngrow=0, grism=grism)
    
    yord, xord = np.indices(orders.shape)
    beams = np.dot(np.ones((orders.shape[0],1), dtype=np.int), xi[5].reshape((1,-1)))
    
    # non_zero = orders > 0
    # 
    # cast = np.int
    # xord, yord, ford, word, sord, bord = np.array(xord[non_zero], dtype=cast), np.array(yord[non_zero], dtype=cast), np.array(orders[non_zero], dtype=np.float64), xi[2][non_zero], xi[3][non_zero], beams[non_zero]
    #         
    ys = orders.shape
    add=1000
    xarr = np.arange(ys[1]+add)+xi[0]
    yarr = np.append(np.sum(orders, axis=0), np.zeros(add))
    
    if grism == 'G141':
        left, right = -195, 84
    else:
        left, right = -195, 84
        
    sum = yarr*0.
    for i in range(1014-left+right):
        sum += nd.shift(yarr, i)
    
    sums = np.zeros((5,ys[1]+add))
    for i in range(5):
        sums[i] = np.append(np.sum(orders*(beams == 2**i), axis=0), np.zeros(add))
    #
    beam_sum = sums*0.
    for i in range(1014-left+right):
        beam_sum += nd.shift(sums, (0, i))
      
    profile = sum[-xi[0]-left:-xi[0]-left+1014]    
    norm = np.median(profile[700:800])
    plt.plot(profile / norm)        
    plt.plot(med_ex / np.median(med_ex[700:800]))
    #plt.plot(med_ex2 / np.median(med_ex2[500:800]))
    
    ### orders
    for i in range(5):
        plt.plot(beam_sum[i, -xi[0]-left:-xi[0]-left+1014] / norm, alpha=0.5, label='BEAM: '+'ABCDE'[i])
        
    
def test_results():
    """
    Run my and Nor's test datasets
    """
    
    roots = ['ibhj02l0q', 'ibhj02leq', 'icdx02vzq', 'icdx03w7q', 'icdx04iwq', 'icdx05jxq', 'icdx05jzq', 'icdx06dbq', 'icdx06ddq', 'icdx07ecq', 'icdx07eeq', 'icdx08enq', 'icdx08epq', 'icdx09j4q']
    
    FORCE=False
    for root in roots:
        if os.path.exists('%s_flt.fits' %(root)) & (not FORCE):
            continue
        #
        #unicorn.prepare.show_MultiAccum_reads('%s_raw.fits' %(root))
        unicorn.prepare.make_IMA_FLT(raw='%s_raw.fits' %(root), pop_reads=[], remove_ima=True)
        
    #
    unicorn.candels.make_asn_files(uniquename=True)
    
    sky_images = {'G141':['zodi_G141_clean.fits', 'excess_lo_G141_clean.fits', 'G141_scattered_light.fits'],
                  'G102':['zodi_G102_clean.fits', 'excess_G102_clean.fits']}
    
    FORCE=False
    
    files=glob.glob('*-G141*asn.fits')
    files=glob.glob('*-G102*asn.fits')
    
    files=files[2:3]
    for file in files:
        #if os.path.exists(file.replace('asn','drz')):
        #    continue
        #
        if 'G102' in file:
            sky = sky_images['G102']
        else:
            sky = sky_images['G141']
        #
        if os.path.exists(file.replace('asn','drz')) & (not FORCE):
            continue
        #
        if not os.path.exists(file.replace('G141', 'F160W').replace('G102', 'F110W')):
            print 'No direct image found: %s' %(file)
            continue
        #
        threedhst.shifts.make_blank_shiftfile(file.replace('G141', 'F140W').replace('G102', 'F110W'), xshift=0, yshift=0, rot=0.0, scale=1.0)
        threedhst.shifts.make_blank_shiftfile(file, xshift=0, yshift=0, rot=0.0, scale=1.0)
        threedhst.prep_flt_files.process_3dhst_pair(file.replace('G141', 'F140W').replace('G102', 'F110W'), file, ALIGN_IMAGE = None, SKIP_GRISM=False, GET_SHIFT=False, SKIP_DIRECT=True, align_geometry='rotate, shift', TWEAKSHIFTS_ONLY=False, adjust_targname=False, sky_images=sky, DIRECT_HIGHER_ORDER=0, GRISM_HIGHER_ORDER=0, final_scale=0.1283, clean_drz=False)
        
    #
    
    
    
          
        
        
              