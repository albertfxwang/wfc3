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
    #new, sub = False, 211
    new, sub = True, 212
    
    ax = fig.add_subplot(sub)

    for root, alpha in zip(['', '-off-x', '-off+x'], [0.4,0.4,0.4]):
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
    
    unicorn.plotting.savefig(fig, 'Vy22_center_Fixed.pdf')
    
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
    
    
    
    
    