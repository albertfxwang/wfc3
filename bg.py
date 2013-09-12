"""
Get spacecraft pointing geometry from Jitter files and compare to grism backgrounds.
"""
import os

import pyfits
import numpy as np
import matplotlib.pyplot as plt

import astropy.time
import astropy.units as u

import datetime

def go():
    import mywfc3.bg
    
    mywfc3.bg.show_orbit_limbangle(asn = ['ib3701050', 'ib3701060'])
    mywfc3.bg.show_orbit_limbangle(asn = ['ib3702050', 'ib3702060'])
    for i in range(28):
        mywfc3.bg.show_orbit_limbangle(asn = ['ib37%02d050' %(i+1), 'ib37%02d060' %(i+1)])
    
    asn_files = glob.glob('ibh*030_*asn.fits')
    for asn in asn_files:
        root = asn.split('_asn')[0]
        mywfc3.bg.show_orbit_limbangle(asn = [root, root.replace('030', '040')])
    
    os.chdir('/Users/brammer/WFC3/Backgrounds/MultiAccum')
    asn_files = glob.glob('*asn.fits')
    for asn in asn_files[7:]:
        root = asn.split('_asn')[0]
        mywfc3.bg.show_orbit_limbangle(asn = [root])
            
    mywfc3.bg.show_orbit_limbangle(asn = ['ibhj44040'])
    
    #mywfc3.bg.show_orbit_limbangle(asn = ['ibhj20030', 'ibhj20040'])

def show_orbit_limbangle(asn = ['ib3701050', 'ib3701060']):
    
    import scipy.ndimage as nd
    
    # os.chdir('/Users/brammer/WFC3/Jitter')
    
    #direct, grism = 'ib3701050', 'ib3701060'
    #direct, grism = 'ibhj20030', 'ibhj20040'
    
    color_list = ['blue','green','red','orange','purple']
    
    jit, colors = [], []
    for i, id in enumerate(asn):
        jit_file = pyfits.open('%s_jit.fits' %(id))
        jit.extend(jit_file[1:])
        colors.extend([color_list[i]]*(len(jit_file)-1))
        
    #jit_direct = pyfits.open('%s_jit.fits' %(direct))
    #jit_grism = pyfits.open('%s_jit.fits' %(grism))
    #jif = pyfits.open(jit.filename().replace('jit','jif'))
    #tstr = time.Time(jif[0].header['TSTRTIME'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
    # jit = jit_direct[1:]
    # jit.extend(jit_grism[1:])
    
    sec_day = u.day.to(u.second)
    
    FLAT_FILE = None
    FLAT_IMAGE = None
    
    FLAT_F140W = pyfits.open(os.path.join(os.getenv('iref'), 'uc721143i_pfl.fits'))[1].data[5:-5, 5:-5]
    
    # colors = ['blue'] * (len(jit_direct)-1)
    # colors.extend(['green'] * (len(jit_grism)-1))
    
    fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_axes((0.05,0.09,0.6,0.28))
    
    ax3 = fig.add_axes((0.05,0.09+0.28,0.6,0.28))
    
    for i, ext in enumerate(jit):
        expname = ext.header['EXPNAME'][:-1]+'q'
        print expname
        spt = pyfits.getheader(expname+'_spt.fits', 0)
        if i == 0:
            targname = spt['TARGNAME']
            ax3.text(0.5, 0.95, targname, ha='center', va='top', transform=ax3.transAxes)
        #
        #### Start/stop times
        pstr = astropy.time.Time(spt['PSTRTIME'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
        pstp = astropy.time.Time(spt['PSTPTIME'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
        if i == 0:
           tstr = pstr
        #
        ax1.text(0.05+0.14*(i % 4), 0.95-0.08*(i/4), expname, ha='left', va='top', color=colors[i], transform=ax1.transAxes, size=9)
        #### Plot curves
        ax1.plot(((pstr-tstr).sec + ext.data['Seconds'])/60., ext.data['LimbAng'], color=colors[i], linewidth=2)
        #
        bright = ext.data['BrightLimb'] == 1
        ax1.plot(((pstr-tstr).sec + ext.data['Seconds'][bright])/60., ext.data['LimbAng'][bright], color=colors[i], linewidth=6, alpha=0.5)
        #
        #day = ext.data['DayNight'] == 0
        #plt.plot(((pstr-tstr).sec + ext.data['Seconds'][bright])/60., ext.data['LimbAng'][bright], color=colors[i], linewidth=12, alpha=0.1)
        #
        #plt.plot((pstr-tstr).sec + ext.data['Seconds'], ext.data['DayNight']*ext.data['TermAng'], color=colors[i], linewidth=3, linestyle=':')
        #
        #plt.plot((pstr-tstr).sec + ext.data['Seconds'], ext.data['DayNight']*20, color=colors[i], linewidth=5, alpha=0.1)
        #
        #### Show background level
        #flt = pyfits.open('../Backgrounds/G141/FLT/%s_flt.fits.gz' %(expname))
        # flt = pyfits.open('%s_flt.fits' %(expname))
        # # if FLAT_FILE != flt[0].header['PFLTFILE']:
        # #     FLAT_IMAGE = pyfits.open(os.path.join(os.getenv('iref'), flt[0].header['PFLTFILE'].split('$')[1]))[1].data[5:-5, 5:-5]
        # #     FLAT_FILE = flt[0].header['PFLTFILE']
        # #
        # flt[1].data /= FLAT_F140W
        # xc, yc, NY = 707, 507, 100
        # subim = flt[1].data[yc-NY:yc+NY, xc-NY:xc+NY]
        # med_background = np.median(subim)
        # ax3.plot(((pstr-tstr).sec + ext.data['Seconds'])/60., ext.data['LimbAng']*0+med_background, color=colors[i], linewidth=2)
        ima = pyfits.open('%s_ima.fits' %(expname))
        time, ramp, reads = get_bg_ramp(ima)
        ax3.plot(((pstr-tstr).sec + time[1:])/60., ramp/np.diff(time), color=colors[i], linewidth=2)
        #
        #img = nd.convolve(flt[1].data, np.ones((2,2))/4.)/med_background
        img = nd.convolve(ima['SCI', 1].data[5:-5, 5:-5], np.ones((2,2))/4.)/np.mean(ramp/np.diff(time))/FLAT_F140W
        if i > -1:
            ax_im = fig.add_axes((0.05+0.30/2*(i % 4)*1.05,0.09+0.28*2+0.01,0.30/2.,0.30))
            ax_im.imshow(img, vmin=0.6, vmax=1/0.6, interpolation='gaussian')
            ax_im.set_xticklabels([]); ax_im.set_yticklabels([])
            trace = np.median(img, axis=0)
            xtr = np.arange(trace.size)
            ax_im.plot(xtr, (trace-1)*2500.+507, color='white', alpha=0.8, linewidth=3)
            ax_im.plot(xtr, (trace-1)*2500.+507, color='black', alpha=0.7, linewidth=1)
            ax_im.set_xlim(0,1014); ax_im.set_ylim(0,1014)
            
    ax1.set_xlabel(r'$\Delta t$ (minutes)')
    ax1.set_ylabel('JIT: LimbAng')
    ax1.set_ylim(0,89)
    xl = ax1.get_xlim()
    ax1.plot([-10,xl[1]], [20,20], linestyle=':')
    ax1.set_xlim(-10, xl[1])
    ax3.set_xlim(-10, xl[1])
    
    ax3.set_ylim(0,3)
    ax3.set_xticklabels([])
    ax3.set_ylabel('backg (electrons/s)')
    ax2 = fig.add_axes((0.66,0.14,0.33,0.4))
    
    
    map = init_map(ax=ax2)
    for i, ext in enumerate(jit):
        draw_map_latlon(map, ext.data['Latitude'], ext.data['Longitude'], color=colors[i], linewidth=2)
        
    tt = np.cast[int](tstr.iso.replace('-',' ').replace(':',' ').replace('.', ' ').split())
    
    #date = datetime.utcnow()
    date = datetime.datetime(tt[0], tt[1], tt[2], tt[3], tt[4], tt[5], tt[6])
    CS=map.nightshade(date, alpha=0.2, color='black')
    ax2.set_title(tstr.iso)
    
    plt.savefig('%s_orbit.png' %(asn[0]))
    plt.close()
    
def init_map(ax=None):
    import mpl_toolkits.basemap as bm
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib.pyplot as plt
    #
    # llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
    # are the lat/lon values of the lower left and upper right corners
    # of the map.
    # lat_ts is the latitude of true scale.
    # resolution = 'c' means use crude resolution coastlines.
    map = Basemap(projection='merc',llcrnrlat=-70,urcrnrlat=70,\
                llcrnrlon=-20, urcrnrlon=380, lat_ts=20, resolution='c', ax=ax)
    #map = Basemap(resolution='c',projection='ortho',lat_0=10.,lon_0=95., ax=ax)
    #
    #m.drawcoastlines()
    map.fillcontinents(color='black',lake_color='aqua', alpha=0.1)
    # draw parallels and meridians.
    #m.drawparallels(np.arange(-90.,91.,30.))
    #m.drawmeridians(np.arange(-180.,181.,60.))
    #map.drawmapboundary(fill_color='aqua')
    map.drawcoastlines(linewidth=0.5)
    map.drawcountries(linewidth=0.5)
    #
    return map
    
def draw_map_latlon(map, lat, lon, *args, **kwargs):
    w = lon < 180
    e = lon >= 180
    xpt, ypt = map(lon[w], lat[w])
    map.plot(xpt, ypt, alpha=0.5, **kwargs)
    xpt, ypt = map(lon[e], lat[e])
    map.plot(xpt, ypt, alpha=0.5, **kwargs)
    xpt, ypt = map(lon[0], lat[0])
    map.scatter(xpt, ypt, alpha=0.5, **kwargs)
    #date = datetime.utcnow()
    #CS=m.nightshade(date, alpha=0.2, color='black')
    
    
def test_ramp(root='ibhj20x7q'):
    import mywfc3.bg
    
    #ima = pyfits.open('%s_ima.fits.gz' %(root))
    ima = pyfits.open('%s_raw.fits' %(root))
    jitfile = glob.glob('%s0[46]0_jit.fits.gz' %(root[:6]))[0]
    jit = pyfits.open(jitfile)
    ext = jit[2]
    
    cube, time, NSAMP = mywfc3.bg.split_multiaccum(ima)
    
    diff = np.diff(cube, axis=0)
    ramp_cps = np.median(diff, axis=1)
    avg_ramp = np.median(ramp_cps, axis=1)
    plt.plot(time[2:], ramp_cps[1:,16:-16:4], alpha=0.1, color='black')
    plt.plot(time[2:], avg_ramp[1:], alpha=0.8, color='red', linewidth=2)
    
    dt = np.diff(time)[3]
    
    #### GOODS-S
    # ok_samp = diff[1:4,:,:]
    
    #ok_samp = diff[-6:,:,:]
    #mi, ma = -3,-1
    #mi, ma = -4, -3
    ok = np.arange(NSAMP-2)[avg_ramp[1:] < 1.2*avg_ramp[1:].min()]
    ok = np.array([np.arange(NSAMP-2)[avg_ramp[1:].argmin()]])
    
    #source = (cube[ma,:,:]-cube[mi,:,:])/(time[ma]-time[mi])*dt
    #source_mean = np.mean(ok_samp, axis=0)
    source = np.sum(diff[ok+1,:,:], axis=0)/(len(ok))

    sky = cube*0.
    for i in range(NSAMP-1):
        sky[i,:,:] = (diff[i,:,:] - source*np.diff(time)[i]/dt)
        
    #sky = cube[:,:,:] - np.median(ok_samp, axis=0)*np.arange(1,NSAMP+1) #*NSAMP #*np.diff(time)[-1]
    for i in range(NSAMP-1):
        ds9.view((sky[i,:,:])/dt)
        
    sky0 = 0.
    for i in range(NSAMP-2):
        sky0 += np.median(sky[i+1,:,:])
    
    sky0 /= (time[-1]-time[1])
    skysum = np.sum(sky[1:,:,:]/(time[-1]-time[1]), axis=0)
    import scipy.ndimage as nd
    skysum_sm = nd.convolve(skysum, np.ones((5,5))/25.)
    cleaned1 = cube[-1,:,:]/time[-1]-skysum
    cleaned0 = cube[-1,:,:]/time[-1]-sky0
    
    #
    FLAT_F140W = pyfits.open(os.path.join(os.getenv('iref'), 'uc721143i_pfl.fits'))[1].data
        
    ramp = np.median(cube, axis=1)
    plt.plot(time[1:], ramp[1:,16:-16:4], alpha=0.1, color='black')
    
    cps = np.diff(ramp[:,16])/np.diff(time)
    plt.scatter(time[1:], cps, color='black')
    
    from scipy import polyfit, polyval
    c = polyfit(time[1:], curve, 3)
    plt.scatter(time[1:], curve, color='black')
    fit = polyval(c, time[1:], color='blue')
    fit0 = polyval(c[-2:], time[1:])
    plt.plot(time[1:], fit0, color='red')
    
def split_multiaccum(ima):
    
    FLAT_F140W = pyfits.open(os.path.join(os.getenv('iref'), 'uc721143i_pfl.fits'))[1].data
    #FLAT_IMAGE = pyfits.open(os.path.join(os.getenv('iref'), ima[0].header['PFLTFILE'].split('$')[1]))[1].data
    
    NSAMP = ima[0].header['NSAMP']
    sh = ima['SCI',1].shape
    cube = np.zeros((NSAMP, sh[0], sh[1]))
    time = np.zeros(NSAMP)
    for i in range(NSAMP):
        #cube[NSAMP-2-i, :, :] = ima['SCI',i+1].data*ima['TIME',i+1].header['PIXVALUE']-ima['SCI',i+2].data*ima['TIME',i+2].header['PIXVALUE']
        if ima[0].header['UNITCORR'] == 'COMPLETE':
            cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data*ima['TIME',i+1].header['PIXVALUE']/FLAT_F140W
        else:
            cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data/FLAT_F140W
        #
        time[NSAMP-1-i] = ima['TIME',i+1].header['PIXVALUE']
    
    return cube, time, NSAMP

### 2D surface http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
import itertools
import numpy as np
import matplotlib.pyplot as plt

def xx():
    """
    Try fitting a 2D polynomial to the sky background
    """
    # Generate Data...
    numdata = 1024
    x = np.random.random(numdata)
    y = np.random.random(numdata)
    z = x**2 + y**2 + 3*x**3 + y + np.random.random(numdata)
    
    yi, xi = np.indices(skysum.shape)
    dq = skysum*0
    dq[5:-5,5:-5] += flt['DQ'].data
    ok = (skysum > 0) & (skysum < 1) & (dq == 0) # & (xi % 2 == 0) & (yi % 2 == 0) & (xi > 10) & (xi < 1014) & (yi > 10) & (yi < 1014) 
    x, y, z = (xi[ok]*1.-512)/512., (yi[ok]*1.-512)/512., skysum[ok]
    
    # Fit a 3rd order, 2d polynomial
    m = polyfit2d(x, y, z, order=33)
    zz = polyval2d((xi*1.-512)/512., (yi*1.-512)/512., m)
    
    fft = np.fft.fft(skysum)
    ifft = np.fft.ifft(fft)
    
    ff = nd.fourier.fourier_gaussian(skysum, 0.1)
    
    # Evaluate it on a grid...
    nx, ny = 20, 20
    xx, yy = np.meshgrid(np.linspace(x.min(), x.max(), nx), 
                         np.linspace(y.min(), y.max(), ny))
    zz = polyval2d(xx, yy, m)

    # Plot
    plt.imshow(zz, extent=(x.min(), y.max(), x.max(), y.min()))
    plt.scatter(x, y, c=z)
    plt.show()

def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

def get_bg_ramp(ima):
    cube, time, NSAMP = split_multiaccum(ima)
    dt = np.diff(time)
    reads = np.diff(cube, axis=0)
    ramp = np.median(reads[:,200:-200,200:-200].reshape((NSAMP-1, -1)), axis=1)
    return time, ramp, reads
    
def test_cal_steps():
    
    FIRST = ['DQICORR', 'ZSIGCORR', 'BLEVCORR', 'ZOFFCORR', 'NLINCORR', 'DARKCORR', 'PHOTCORR']
    SECOND = ['UNITCORR', 'CRCORR', 'FLATCORR']

    FIRST = ['DQICORR', 'ZSIGCORR', 'BLEVCORR', 'ZOFFCORR', 'NLINCORR', 'DARKCORR', 'PHOTCORR', 'UNITCORR', 'FLATCORR']
    SECOND = ['CRCORR']
    
    file='ibhj44mrq_raw.fits'
    mywfc3.bg.run_wf3ir(file=file, PERFORM=FIRST, OMIT=SECOND)
    
    ima = pyfits.open(file.replace('raw', 'ima'), mode='update')
    cube, time, NSAMP = mywfc3.bg.split_multiaccum(ima)
    dt = np.diff(time)
    reads = np.diff(cube, axis=0)
    ramp = np.median(reads[:,200:-200,200:-200].reshape((NSAMP-1, -1)), axis=1)
    ok = dt > 50

    cut_ramp = np.median(reads[:,5:-5,5:-5], axis=1)
    
    ramp_counts = np.median(reads[:,5:-5,:], axis=1)
    avg_ramp = np.median(ramp_counts[:,5:-5], axis=1)
    plt.plot(time[2:], ramp_counts[1:,16:-16:4], alpha=0.1, color='black')
    plt.plot(time[2:], avg_ramp[1:], alpha=0.8, color='red', linewidth=2)
    
    min_bg_idx = avg_ramp[1:].argmin()+1
    baseline_object_cps = reads[min_bg_idx,:,:]/dt[min_bg_idx]
    baseline_object_dq = ima['DQ', NSAMP-min_bg_idx-1].data
    
    sky = cube*0.
    for i in range(NSAMP-1):
        sky[i+1,:,:] = (reads[i,:,:] - baseline_object_cps*dt[i])
    
    avg_excess = np.median(sky[:,5:-5,5:-5].reshape((NSAMP, -1)), axis=1)
    sum_excess = np.cumsum(avg_excess)
    for i in range(NSAMP):
        ima['SCI', i+1].data -= sum_excess[-(i+1)]
    #
    ima.flush()
    #
    # for i in range(NSAMP-1):
    #     ds9.view((sky[i+1,:,:])/dt)
    
    mywfc3.bg.run_wf3ir(file='ibhj44mrq_ima.fits', PERFORM=SECOND, OMIT=FIRST)
    
def run_wf3ir(file='ibhj44mrq_raw.fits', 
              OMIT = ['UNITCORR', 'CRCORR', 'FLATCORR'], 
              PERFORM = ['DQICORR', 'ZSIGCORR', 'BLEVCORR', 'ZOFFCORR', 'NLINCORR', 'DARKCORR', 'PHOTCORR'], 
              verbose=True, clean=True, 
              ):
    """
    Perform steps before up-the-ramp sampling on a raw image to improve time-variable
    background subtraction of extra earthglow
    """
    os.chdir("/Users/brammer/WFC3/Backgrounds/CalWFC3")
    from wfc3tools import calwf3
    
    im = pyfits.open(file, mode='update')
    for key in PERFORM:
        im[0].header.update(key, 'PERFORM')
    
    for key in OMIT:
        im[0].header.update(key, 'OMIT')
    
    im.flush()
    
    flt = file.replace('ima', 'flt').replace('raw', 'flt')
    
    if clean:
        if 'raw' in file:
            for ext in ['ima', 'flt']:
                if os.path.exists(file.replace('raw', ext)):
                    os.remove(file.replace('raw', ext))
        else:
            for ext in ['ima', 'flt']:
                if os.path.exists(file.replace('ima', 'ima_%s' %ext)):
                    os.remove(file.replace('ima', 'ima_%s' %ext))
            
    calwf3.calwf3(file, verbose=verbose)
    
    
    
    
    
    