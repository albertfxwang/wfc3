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
    
    ima = pyfits.open(files[ix])
    grism = ima[0].header['FILTER']
    imaging_flats = {'G102':'uc72113oi_pfl.fits', 'G141':'uc721143i_pfl.fits'} ## F105/F140W
    imaging_flats = {'G102':'uc72113ni_pfl.fits', 'G141':'uc721143i_pfl.fits'} ## F098M/F140W
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
    pyfits.writeto('%s_earth.fits' %(ima.filename().split('_ima')[0]), data=sum, header=ima['SCI',1].header, clobber=True)
    
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

    sky = pyfits.open(sky_file[grism])[0].data/PAM/flat#*grism_flat
    if 'set' in sky_file[grism]: sky*=flat
    
    sky_med = np.median(sky[100:900,400:600])
    ds9.frame(i+3)
    sky[sky == 0] = sum[sky == 0]/sum_med
    ds9.view(sky)
    ds9.set('scale limits %f %f' %(sky_med*0.7, sky_med*1.3))
    
    ds9.set('lock frame image'); ds9.set('lock colorbar')
    ds9.set('cmap hsv'); ds9.set('width 1335'); ds9.set('height 890'); ds9.set('cmap value 4.86047 0.507116')
    ds9.set('tile yes'); ds9.set('zoom to fit')
    
    #ds9.set('saveimage gif test.gif')
    ds9.set('saveimage tiff %s.tiff none' %(ima.filename().split('_ima')[0]))
    
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
    
    