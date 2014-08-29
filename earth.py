import pyfits

def prep_test():
    import mywfc3.bg
    import glob
    
    files=glob.glob('*raw.fits')
    files=glob.glob('ich211*raw.fits')
    for file in files:
        mywfc3.bg.run_wf3ir(file=file, OMIT = ['CRCORR'], PERFORM = ['DQICORR', 'ZSIGCORR', 'BLEVCORR', 'ZOFFCORR', 'NLINCORR', 'DARKCORR', 'PHOTCORR', 'UNITCORR', 'FLATCORR'], verbose=True, clean=True)

    files=glob.glob('*ima.fits')
    #files=glob.glob('ich211*ima.fits')
    ix = 0
    ix+=1
    
    ima = pyfits.open(files[ix])
    grism = ima[0].header['FILTER']
    imaging_flats = {'G102':'uc72113oi_pfl.fits', 'G141':'uc721143i_pfl.fits'} ## F105/F140W
    #imaging_flats = {'G102':'uc72113ni_pfl.fits', 'G141':'uc721143i_pfl.fits'} ## F098M/F140W
    #imaging_flats = {'G102':'uc72113ni_pfl.fits', 'G141':'uc721145i_pfl.fits'} ## F098M/F125W

    grism_flats = {'G102':'u4m1335li_pfl.fits', 'G141':'u4m1335mi_pfl.fits'} ## F098M/F140W
    
    image_flat = pyfits.open(os.getenv('iref')+'/'+imaging_flats[grism])[1].data[5:-5,5:-5]
    grism_flat = pyfits.open(os.getenv('iref')+'/'+grism_flats[grism])[1].data[5:-5,5:-5]
    
    flat = image_flat/grism_flat
    #flat = 1
    
    cube, time, NSAMP = mywfc3.bg.split_multiaccum(ima, scale_flat=False)
    
    diff = np.diff(cube, axis=0)
    dt = np.diff(time)
    ramp_cps = np.median(diff, axis=1)
    avg_ramp = np.median(ramp_cps, axis=1)
    
    PAM = pyfits.open(os.getenv('iref')+'/ir_wfc3_map.fits')[1].data
    PAM = 1
    #ds9 = threedhst.dq.myDS9()
    for i in range(len(dt)-1):
        ds9.frame(i+1)
        corr = diff[i+1,5:-5,5:-5]/PAM/flat
        med = np.median(corr[100:900,400:600])
        ds9.view(corr)
        ds9.set('scale limits %f %f' %(med*0.7, med*1.3))
        
    sum = diff[1:,5:-5,5:-5].sum(axis=0)/PAM/flat
    sum_med = np.median(sum[100:900,400:600])
    ds9.frame(i+2)
    ds9.view(sum)
    ds9.set('scale limits %f %f' %(sum_med*0.7, sum_med*1.3))
    
    ### Correct for read pattern 
    medx = np.dot(np.median(sum, axis=0).reshape((1014,1)), np.ones((1,1014))).T
    medy = np.dot(np.median(sum/medx, axis=1).reshape((1014,1)), np.ones((1,1014)))
    # normed = sum / medx
    # for i in range(1014)[::10]:
    #     plt.plot(normed[:,i], alpha=0.01, color='blue')
    sum /= medy
    flt = pyfits.open(files[ix].replace('_ima','_flt')+'.gz')
    flt[1].data = sum
    flt.writeto('%s_earth.fits' %(ima.filename().split('_ima')[0]), clobber=True)
    
    #zodi = pyfits.open('/Users/brammer/3DHST/Spectra/Work/CONF/zodi_G102_clean.fits')[0].data
    zodi = pyfits.open('/Users/brammer/3DHST/Spectra/Work/CONF/zodi_%s_clean.fits' %(grism))[0].data
    #flt[1].data /= zodi
    #flt.writeto('%s_x.fits' %(ima.filename().split('_ima')[0]), clobber=True)
    
    sum_med = np.median(sum[100:900,400:600])
    ds9.view(sum)
    ds9.set('scale limits %f %f' %(sum_med*0.7, sum_med*1.3))
    
    fp = open('label.reg','w')
    fp.writelines(['image\n', '# text(507.68945,886.64258) color=white font="helvetica 12 bold roman" text={%s: %s, %d electrons}\n' %(grism, ima.filename().split('_ima')[0], sum_med)])
    fp.close()
    ds9.set('regions delete all')
    ds9.set('regions file label.reg')
    
    sky_file = {'G102': os.getenv('THREEDHST')+'/CONF/WFC3.IR.G102.sky.V1.0.fits', 
                'G141': os.getenv('THREEDHST')+'/CONF/WFC3.IR.G141.sky.V1.0.fits'}
    
    #sky_file['G141'] = '/Users/brammer/3DHST/Spectra/Work/CONF/sky.G141.set002.fits'

    #sky = pyfits.open(sky_file[grism])[0].data/PAM/flat#*grism_flat
    #if 'set' in sky_file[grism]: sky*=flat
    sky = zodi
    
    sky_med = np.median(sky[100:900,400:600])
    ds9.frame(i+3)
    sky[sky == 0] = sum[sky == 0]/sum_med
    ds9.view(sky)
    ds9.set('scale limits %f %f' %(sky_med*0.7, sky_med*1.3))
    
    ds9.set('lock frame image'); ds9.set('lock colorbar')
    ds9.set('cmap hsv'); ds9.set('width 1100'); ds9.set('height 725'); ds9.set('cmap value 4.86047 0.507116')
    ds9.set('tile yes'); ds9.set('zoom to fit')
    
    #ds9.set('saveimage gif test.gif')
    #ds9.set('saveimage tiff %s_%s.tiff none' %(ima.filename().split('_ima')[0], grism))
    
    import scipy.ndimage as nd
    
    fig = unicorn.plotting.plot_init(aspect=2./3, xs=10, left=0.01, right=0.01, top=0.01, bottom=0.01, square=True, hspace=0., wspace=0.)
    
    vm = [0.95, 1.05]
    
    for i in range(4):
        corr = diff[i+1,5:-5,5:-5]/PAM/flat
        med = np.median(corr[100:900,400:600])
        ax = fig.add_subplot(231+i)
        ax.imshow(nd.gaussian_filter(corr, 1), interpolation='Nearest', vmin=med*0.95, vmax=med*1.05)
        #ds9.view(corr)
        #ds9.set('scale limits %f %f' %(med*0.9, med*1.3))
    
    ax = fig.add_subplot(235)
    med = np.median(sum[100:900,400:600])
    ax.imshow(nd.gaussian_filter(sum, 1), interpolation='Nearest', vmin=med*0.95, vmax=med*1.05)
    ax.text(0.5, 0.95, '%s %s' %(ima.filename().split('_ima')[0], grism), ha='center', va='top', transform=ax.transAxes, backgroundcolor='white')
    ax.text(0.5, 0.05, '%d e-' %(sum_med), ha='center', va='bottom', transform=ax.transAxes, backgroundcolor='white')
    
    ax = fig.add_subplot(236)
    med = np.median(sky[100:900,400:600])
    ax.imshow(nd.gaussian_filter(sky, 1), interpolation='Nearest', vmin=med*0.95, vmax=med*1.05)
    
    for ax in fig.axes:
        ax.set_xticks([0,1014]); ax.set_yticks([0,1014])
        ax.set_xticklabels([]); ax.set_yticklabels([])
        
    unicorn.plotting.savefig(fig, '%s_%s.png' %(ima.filename().split('_ima')[0], grism))
    
    ######## Histograms
    import scipy.ndimage as nd
    kernel = np.ones((10,10))
    local_mean_sum = nd.convolve(sum, kernel)/kernel.sum()
    resid_sum = (sum-local_mean_sum)/local_mean_sum

    local_mean_sky = nd.convolve(sky, kernel)/kernel.sum()
    resid_sky = (sky-local_mean_sky)/local_mean_sky
    
    plt.hist(resid_sum.flatten(), range=(-0.05,0.05), bins=100, alpha=0.3, label='DARK-EARTH')
    plt.hist(resid_sky.flatten(), range=(-0.05,0.05), bins=100, alpha=0.3, label='On-sky')
    plt.legend()
    plt.xlim(-0.05,0.05)
    
    ### Use scripts in bg_ISR to identify the strong (OH?) line in G141
    im = pyfits.open('../RAW/ibhj03xvq_flt.fits', mode='update')
    new = pyfits.open('ich213iuq_earth.fits')
    im[1].data = new[0].data
    im.flush()
    
def extract_blob(xc=231, yc=417, grism='G141', NY=40, file='ich219i8q_earth.fits', ext=1, show=True):
    """
    extract spectrum of a blob
    """
    import unicorn
    import os
    import pyfits
    import numpy as np
    import matplotlib.pyplot as plt
    
    orders, xi = unicorn.reduce.grism_model(xc, yc, lam_spec=None, flux_spec=None, BEAMS=['A'], grow_factor=1, growx=1, growy=1, pad=0, ngrow=1, grism=grism)
    
    yord, xord = np.indices(orders.shape)
    beams = np.dot(np.ones((orders.shape[0],1), dtype=np.int), xi[5].reshape((1,-1)))
    
    non_zero = orders > 0
    
    cast = np.int
    xord, yord, ford, word, sord, bord = np.array(xord[non_zero], dtype=cast), np.array(yord[non_zero], dtype=cast), np.array(orders[non_zero], dtype=np.float64), xi[2][non_zero], xi[3][non_zero], beams[non_zero]
            
    ys = orders.shape
    xord += xi[0]
    yord -= (ys[0]-1)/2
    #        
    xxi = np.int(np.round(xc))+xord
    use = (xxi >= 0) & (xxi < 1014)
    ### First order only
    #use = use & (xord > 10*self.grow_factor) & (xord < 213*self.grow_factor) 
    use_order = use 
    cut = np.zeros(1014)
    cut[xxi[use_order]] = word[use_order] 
    object_wave = cut.copy() #np.dot(np.ones((self.sh[0],1)),
    
    wlim = unicorn.reduce.grism_wlimit[grism]
    
    wavelength_region = (object_wave >= wlim[0]*1) & (object_wave <= wlim[1]*1.)
    
    data = pyfits.open(file)[ext].data
    dq = pyfits.open(file)['DQ'].data
        
    #### dydx offset for grism trace
    xarr = np.arange(1014)
    xmin = xarr[wavelength_region].min()
    xmax = xarr[wavelength_region].max()
    
    xoff_arr = np.arange(len(xi[4]), dtype=np.double)+xi[0]+xc
    xoff_int = np.arange(xmin, xmax+0.1, dtype=np.double)
    yoff_int = np.interp(xoff_int, xoff_arr, xi[4])
    
    wave = object_wave[wavelength_region]
    y0 = int(yc+np.interp(wlim[2], wave, yoff_int))
    
    twod = data[y0-NY:y0+NY, wavelength_region]
    twod_dq = dq[y0-NY:y0+NY, wavelength_region] == 0
    twod[twod_dq == 0] = 0
    
    oned = np.sum(twod[NY-10:NY+10,:], axis=0)
    oned_dq = np.sum(twod_dq[NY-10:NY+10,:], axis=0)
    oned = oned / oned_dq
    
    bg = np.sum(twod[2*NY-20:,:], axis=0)
    bg_dq = np.sum(twod_dq[2*NY-20:,:], axis=0)
    bg = bg / bg_dq
        
    unicorn.reduce.set_grism_config(grism=grism)
    sens = pyfits.open(os.getenv('THREEDHST')+'/CONF/'+unicorn.reduce.conf['SENSITIVITY_A'])[1].data
    sens_interp = np.interp(wave, sens['WAVELENGTH'], sens['SENSITIVITY'])
    sens_interp /= sens_interp.max()
    
    full_spec = (bg/oned-1)/sens_interp
    if show:
        plt.plot(wave, full_spec, alpha=0.2, linewidth=2, color='black')
    
    return wave, full_spec, twod
    
def show_all():
    
    files=glob.glob('*earth.fits')
    xc, yc = 230, 417
    #xc, yc = 468, 374
    #xc, yc = 127, 964
    
    NY=40
    
    file = files[20]
    im = pyfits.open(file)
    grism = im[0].header['FILTER']
    wave, full_spec, twod = mywfc3.earth.extract_blob(file=file, grism=grism, show=False, xc=xc, yc=yc, NY=NY)
        
    bad = ['ich203g9q', 'ich204ibq', 'ich204icq', 'ich207wiq', 'ich210teq', 'ich210tfq', 'ich212hiq', 'ich212hkq', 'ich213isq', 'ich217sbq', 'ich219i6q']

    # im = pyfits.open(files[0])[1].data
    # thumb = im[yc-NY:yc+NY, xc-NY:xc+NY]
    # #thumb /= np.median(thumb[:,0:10])
    # 
    # kernel = np.sum(thumb[NY-10:NY+10,:], axis=0)
    # kernel = kernel - np.sum(thumb[-20:,:], axis=0)
    # kernel[kernel < 0] = 0
    # kernel /= np.sum(kernel)
    # ywave = wave*0
    # ywave[59] = 1
    # yy = nd.convolve(ywave, kernel)
    
    sum = full_spec*0.
    
    count = 0
    for file in files[20:]: 
        if file.split('_')[0] in bad:
            continue
        #
        count += 1
        #
        im = pyfits.open(file)
        grism = im[0].header['FILTER']
        wave, full_spec, twod = mywfc3.earth.extract_blob(file=file, grism=grism, show=False, xc=xc, yc=yc, NY=NY)
        #
        ds9.view(twod/np.median(twod))
        plt.plot(wave, full_spec, alpha=0.2, linewidth=2, color='black')
        sum += full_spec
    
    plt.plot(wave, sum/count, color='red', alpha=0.8, linewidth=2)
    
    np.savetxt('dark_earth_G141.dat', np.array([wave, sum/count]).T, fmt='%13.5e')
    
    plt.ylim(-0.005, 0.1)
    
        
    