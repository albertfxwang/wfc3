"""
Flux calibration
"""
import numpy as np
import matplotlib.pyplot as plt

def check_GD153():
    import mywfc3.flux
    
    os.chdir('/Users/brammer/WFC3/Calibration/Cycle20/13092_Flux/Reduce')
    
    asn = threedhst.utils.ASNFile('Vy22-blue-F098M_asn.fits')
    
    root=''
    images = {'F098M':'ic461aghq_flt.fits.gz', 'F105W':'ic461agiq_flt.fits.gz', 'G102':'ic461agjq_flt.fits.gz', 'F140W':'ic461ah4q_flt.fits.gz', 'F160W':'ic461ah5q_flt.fits.gz', 'G141':'ic461ah6q_flt.fits.gz'}

    root='-off+x'
    images = {'F098M':'ic461agpq_flt.fits.gz', 'F105W':'ic461agqq_flt.fits.gz', 'G102':'ic461agrq_flt.fits.gz', 'F140W':'ic461ahbq_flt.fits.gz', 'F160W':'ic461ahcq_flt.fits.gz', 'G141':'ic461ahdq_flt.fits.gz'}
    # 
    # root='-off-x'
    # images = {'F098M':'ic5v02aiq_flt.fits.gz', 'F105W':'ic5v02ajq_flt.fits.gz', 'G102':'ic5v02akq_flt.fits.gz', 'F140W':'ic5v41b6q_flt.fits.gz', 'F160W':'ic5v41b7q_flt.fits.gz', 'G141':'ic5v41b8q_flt.fits.gz'}

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
        asn.product = 'GD153%s-%s-%s' %(root, band, filter)
        asn.exposures = [images[filter].split('_flt')[0]]
        asn.write(asn.product+'_asn.fits')
        im = pyfits.open('../RAW/'+images[filter])
        if filter in flat.keys():
            im[1].data /= flat[key]
            sky = pyfits.open('/Users/brammer/3DHST/Spectra/Work/CONF/sky.G141.set002.fits ')[0].data
            #sky /= flat[key]
            a = np.median(im[1].data/sky)
            im[1].data -= a*sky
        else:
            im[1].data -= np.median(im[1].data)
        #
        #if not os.path.exists(images[filter][:-3]):
        im.writeto(images[filter].split('.gz')[0], clobber=True)
    
    files=glob.glob('GD*%s*asn.fits' %(root))
    for file in files:
        unicorn.reduce.interlace_combine(file.split('_asn')[0], growx=1, growy=1, pad=60, NGROW=0, view=False)
        
    model = unicorn.reduce.GrismModel('GD153-blue', direct='F105W', grism='G102', growx=1, growy=1, grow_factor=1)
    model.twod_spectrum(20, miny=-50)
    model.show_2d(savePNG=True)
    
    model = unicorn.reduce.GrismModel('GD153-red', growx=1, growy=1, grow_factor=1)
    model.twod_spectrum(12, miny=-50)
    model.show_2d(savePNG=True)
    
    import scipy.ndimage as nd

    twod = unicorn.reduce.Interlace2D('GD153-red_00012.2D.fits')
    wave, flux = twod.optimal_extract(twod.im['SCI'].data-twod.im['CONTAM'].data)
    sens = twod.im['SENS'].data    
    plt.plot(wave, flux/sens, color='black', alpha=0.5, linewidth=3)

    plt.plot(wave, flux/nd.shift(sens, -0.5), color='purple', alpha=0.5)
    #plt.plot(wave, flux/nd.shift(sens, -1), color='green', alpha=0.5)
    # plt.plot(wave, flux/nd.shift(sens, -1.), color='red', alpha=0.5)
    
    sp = pyfits.open('/grp/hst/cdbs/calspec/gd153_mod_008.fits')[1].data
    plt.plot(sp['WAVELENGTH'], sp['FLUX']/1.e-17*0.92)
    
    twod.compute_model(lam_spec=np.cast[float](sp['WAVELENGTH']), flux_spec=np.cast[float](sp['FLUX'])/1.e-17*0.92/twod.total_flux)
    twod.model[twod.im['SCI'].data == 0] = 0
    w2, f2 = twod.optimal_extract(twod.model)
    plt.plot(w2, f2/sens, color='green')

    bl = unicorn.reduce.Interlace2D('GD153-blue_00020.2D.fits')
    blwave, blflux = bl.optimal_extract(bl.im['SCI'].data-bl.im['CONTAM'].data)
    blsens = bl.im['SENS'].data    
    plt.plot(blwave, blflux/blsens, color='black', alpha=0.5, linewidth=3)
    plt.plot(blwave, blflux/nd.shift(blsens, -0.5), color='purple', alpha=0.5)

    bl.compute_model(lam_spec=np.cast[float](sp['WAVELENGTH']), flux_spec=np.cast[float](sp['FLUX'])/1.e-17*0.92/bl.total_flux)
    bl.model[bl.im['SCI'].data == 0] = 0
    w2, f2 = bl.optimal_extract(bl.model)
    plt.plot(w2, f2/blsens, color='green')

    
    plt.xlim(0.7e4,1.7e4)
    plt.ylim(10,500)
    plt.semilogy()
    
    model = unicorn.reduce.GrismModel('GD153-off+x-blue', direct='F105W', grism='G102', growx=1, growy=1, grow_factor=1)
    model.twod_spectrum(17, miny=-50)
    model.show_2d(savePNG=True)
    
    twod = unicorn.reduce.Interlace2D('GD153-off+x-blue_00017.2D.fits')
    w, f = twod.optimal_extract(twod.im['SCI'].data-twod.im['CONTAM'].data)
    s = twod.im['SENS'].data    
    plt.plot(w, f/s, color='red', alpha=0.5, linewidth=3)
    
    model = unicorn.reduce.GrismModel('GD153-off+x-red', direct='F140W', grism='G141', growx=1, growy=1, grow_factor=1)
    model.twod_spectrum(18, miny=-50)
    model.show_2d(savePNG=True)
    
    twod = unicorn.reduce.Interlace2D('GD153-off+x-red_00018.2D.fits')
    w, f = twod.optimal_extract(twod.im['SCI'].data-twod.im['CONTAM'].data)
    s = twod.im['SENS'].data    
    plt.plot(w, f/s, color='red', alpha=0.5, linewidth=3)
    
    