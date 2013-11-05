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

from mywfc3.utils import gzfile
import unicorn

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
    asn_files = glob.glob('*jit.fits*')
    for asn in asn_files:
        root = asn.split('_jit')[0]
        mywfc3.bg.show_orbit_limbangle(asn = [root])
        
    mywfc3.bg.show_orbit_limbangle(asn = ['ibhj44040'])
    
    ### Overall programs
    os.chdir('/Users/brammer/WFC3/GrismPrograms/CANDELS-SNe/RAW/')
    
    SKIP=True
    
    asn_files = glob.glob('*asn.fits')
    for asn in asn_files:
        root = asn.split('_asn')[0]
        png = glob.glob(root+'*orbit.*')
        if (len(png) > 0) & SKIP:
            continue
        try:
            mywfc3.bg.show_orbit_limbangle(asn = [root])
        except:
            fp = open(root+'_orbit.failed','w')
            fp.write('---')
            fp.close()
            
    #mywfc3.bg.show_orbit_limbangle(asn = ['ibhj20030', 'ibhj20040'])

def show_orbit_limbangle(asn = ['ib3701050', 'ib3701060']):
    
    import scipy.ndimage as nd
    
    # os.chdir('/Users/brammer/WFC3/Jitter')
    
    #direct, grism = 'ib3701050', 'ib3701060'
    #direct, grism = 'ibhj20030', 'ibhj20040'
    
    color_list = ['blue','green','red','orange','purple']
    
    jit, colors = [], []
    for i, id in enumerate(asn):
        jit_file = pyfits.open(gzfile('%s_jit.fits' %(id)))
        jit.extend(jit_file[1:])
        colors.extend([color_list[i]]*(len(jit_file)-1))
    
    try:
        jif = pyfits.open(gzfile('%s_jif.fits' %(asn[0])))[1:]
    except:
        jif = None
        
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
    FLAT_F105W = pyfits.open(os.path.join(os.getenv('iref'), 'uc72113oi_pfl.fits'))[1].data[5:-5, 5:-5]
    
    # colors = ['blue'] * (len(jit_direct)-1)
    # colors.extend(['green'] * (len(jit_grism)-1))
    
    #fig = plt.figure(figsize=(12,6))
    fig = unicorn.plotting.plot_init(xs=12, aspect=0.5, NO_GUI=True)
    ax1 = fig.add_axes((0.05,0.09,0.6,0.28))
    
    ax3 = fig.add_axes((0.05,0.09+0.28,0.6,0.28))
    
    for i, ext in enumerate(jit):
        expname = ext.header['EXPNAME'][:-1]+'q'
        print expname
        spt = pyfits.getheader(gzfile(expname+'_spt.fits'), 0)
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
        ok = ext.data['TermAng'] < 1000
        ax1.plot((((pstr-tstr).sec + ext.data['Seconds'])/60.)[ok], ext.data['TermAng'][ok], color='green', alpha=0.4, linewidth=2)
        ax1.plot(((pstr-tstr).sec + ext.data['Seconds'])/60., ext.data['LOS_Zenith'], color='orange', alpha=0.4, linewidth=2)
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
        ima = pyfits.open(gzfile('%s_raw.fits' %(expname)))
        FILTER = ima[0].header['FILTER']
        time, ramp, reads = get_bg_ramp(ima)
        dt = np.diff(time)
        ok = dt > 24
        ax3.plot(((pstr-tstr).sec + time[1:][ok])/60., (ramp/dt*2.5)[ok], color=colors[i], linewidth=2)
        #
        ##
        ## JIF Shadow
        if jif is not None:
            shadow_in = astropy.time.Time(jif[i].header['SHADOENT'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
            shadow_out = astropy.time.Time(jif[i].header['SHADOEXT'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
            dt_in = (shadow_in-pstr).sec
            dt_out = (shadow_out-pstr).sec
            print 'JIF: in=%.3f out=%.3f' %(dt_in, dt_out)
            if dt_in > 0:
                yi = np.interp(dt_in, (time[1:][ok]), (ramp/dt*2.5)[ok])
                ax3.scatter(((pstr-tstr).sec + dt_in)/60., yi, marker='o', s=40, color='green', alpha=0.8)
            #
            if (dt_out > 0) & (dt_out < time.max()):
                yi = np.interp(dt_out, (time[1:][ok]), (ramp/dt*2.5)[ok])
                ax3.scatter(((pstr-tstr).sec + dt_out)/60., yi, marker='o', s=40, color='red', alpha=0.8)
                
        #### Make log
        idx = np.arange(len(time)-1)[ok]
        fp = open('%s_%s_orbit.dat' %(ext.header['EXPNAME'], FILTER), 'w')
        dtypes = {'E':'%.3f', 'D':'%.3f', 'I':'%d'}
        line = '# name read  bg minbg'
        for c in ext.data.columns:
            line += ' %s' %(c.name)
        fp.write(line+'\n')
        
        for il in idx:
            line = '%s_%02d %2d %.3f %.3f' %(ext.header['EXPNAME'], il, il, ramp[il]/dt[il]*2.5, np.min((ramp/dt*2.5)[ok]))
            for c in ext.data.columns:
                val = np.interp(time[1:][il]-dt[il]/2., ext.data['Seconds'], ext.data[c.name])
                line += ' ' + dtypes[c.format] %(val) # %(str(np.cast[dtypes[c.format]](val)))
            #
            fp.write(line+'\n')
        fp.close()
        
        #img = nd.convolve(flt[1].data, np.ones((2,2))/4.)/med_background
        last_read = (ima['SCI', 1].data-ima['SCI',len(time)].data)[5:-5, 5:-5]*2.5
        FLAT = FLAT_F140W
        if FILTER in ['G102', 'F105W']:
            FLAT = FLAT_F105W
            
        last_read = nd.minimum_filter(last_read/FLAT, size=2)
        img = nd.convolve(last_read, np.ones((2,2))/4.)/time[-1]
        img /= np.median(img[400:600,400:600])
        
        if i > -1:
            ax_im = fig.add_axes((0.05+0.30/2*i*1.05,0.09+0.28*2+0.01,0.30/2.,0.30))
            ax_im.imshow(img, vmin=0.65, vmax=1/0.57, interpolation='gaussian')
            ax_im.set_xticklabels([]); ax_im.set_yticklabels([])
            trace = np.median(img, axis=0)
            xtr = np.arange(trace.size)
            ax_im.plot(xtr, (trace-1)*2500.+507, color='white', alpha=0.8, linewidth=3)
            ax_im.plot(xtr, (trace-1)*2500.+507, color='black', alpha=0.7, linewidth=1)
            ax_im.set_xlim(0,1014); ax_im.set_ylim(0,1014)
            
    ax1.set_xlabel(r'$\Delta t$ (minutes)')
    ax1.set_ylabel('JIT: LimbAng')
    ax1.set_ylim(0,99)
    xl = ax1.get_xlim()
    ax1.plot([-10,xl[1]], [20,20], linestyle=':')
    ax1.set_xlim(-10, xl[1])
    ax3.set_xlim(-10, xl[1])
    #ax3.legend(loc='upper left', prop={'size':8})
    
    ax3.set_ylim(0,3.8)
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
    
    #plt.savefig('%s_%s_orbit.png' %(asn[0], FILTER))
    unicorn.plotting.savefig(fig, '%s_%s_orbit.png' %(asn[0], FILTER))
    #plt.close()
    
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
    #os.chdir("/Users/brammer/WFC3/Backgrounds/CalWFC3")
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
    
    
def background_emission_line():
    """
    In some high background G141 observations, it looks like the dispersed "blobs" become point sources indicative
    of a stong emission-line component of the background.  That feature is at the far blue edge of the G141 spectrum 
    and therefore difficult to determine the wavelength for.  
    
    It turns out that the feature can also be seen in G102, specifically exposure 'ibl004lwq' from Koekemoer's AGN 
    program.  Figure out what the wavelength is.
    
    He I, a strong line in the sunlit aurora
    
    """
    
    os.chdir('/Users/brammer/WFC3/GrismPrograms/Koekemoer/BackgroundG102')
    
    asn = threedhst.utils.ASNFile('../RAW/ibl004010_asn.fits')
    
    asn.exposures = ['ibl004lwq']
    asn.product = 'BLOB-G102'
    asn.write('BLOB-G102_asn.fits')

    asn.exposures = ['ibl004l3q']
    asn.product = 'BLOB-F104W'
    asn.write('BLOB-F140W_asn.fits')
    
    threedhst.process_grism.fresh_flt_files('BLOB-G102_asn.fits')
    
    #### Subtract background and add a source for the blob extraction
    FLAT_F105W = pyfits.open(os.path.join(os.getenv('iref'), 'uc72113oi_pfl.fits'))[1].data[5:-5, 5:-5]
    blobs = [(230, 418), (467, 377)]
    source = FLAT_F105W*0.
    for blob in blobs:
        yi, xi = np.indices(FLAT_F105W.shape)
        r = np.sqrt((xi-blob[0])**2+(yi-blob[1])**2)
        source[r < 20] += 1./FLAT_F105W[r < 20]-1
        #source[r > 20] = 0
    
    threedhst.process_grism.fresh_flt_files('BLOB-F140W_asn.fits')
    im = pyfits.open('ibl004l3q_flt.fits', mode='update')
    im[1].data += source*50-np.median(im[1].data)
    im.flush()

    unicorn.reduce.interlace_combine('BLOB-F140W', growx=1, growy=1)
    threedhst.shifts.make_blank_shiftfile('BLOB-F140W_asn.fits')
    threedhst.prep_flt_files.startMultidrizzle('BLOB-F140W_asn.fits', final_scale=0.128254, pixfrac=1)
    
    unicorn.reduce.interlace_combine('BLOB-G102', growx=1, growy=1)
    
    model = unicorn.reduce.GrismModel(root='BLOB', grow_factor=1, growx=1, growy=1, MAG_LIMIT=21, use_segm=False, grism='G102')
    ### BLob
    model.twod_spectrum(108, miny=-100)
    model.twod_spectrum(119, miny=-100)
    ### Star
    model.twod_spectrum(74, miny=-100)
    
    twod = unicorn.reduce.Interlace2D('BLOB_00108.2D.fits')
    spec = twod.im['SCI'].data
    spec = np.median(spec)/spec-1
    lam, spec1d = twod.optimal_extract(spec)
    
    kernel = np.sum(twod.thumb**2, axis=0)
    kernel -= np.mean(kernel[0:30])
    kernel /= kernel.sum()

    twod = unicorn.reduce.Interlace2D('BLOB_00119.2D.fits')
    spec = twod.im['SCI'].data
    spec = np.median(spec)/spec-1
    lam_B, spec1d_B = twod.optimal_extract(spec)
    
    cc = nd.correlate1d(spec1d, kernel)
    cc_B = nd.correlate1d(spec1d_B, kernel)
    wavelen = lam[np.argmax(cc)] ## 1.083 um = He I !
    peak = np.arange(np.argmax(cc)-4,np.argmax(cc)+8)
    cc_peak = np.trapz(lam[peak]*cc[peak], lam[peak]) / np.trapz(cc[peak], lam[peak])
    
    plt.plot(lam[0:kernel.shape[0]], kernel/kernel.max()*2, label='Inverted blob (F105W)', linestyle='steps')
    plt.plot(lam, spec1d, label='Blob spectrum (230, 418)', linestyle='steps')
    plt.plot(lam, cc, label=r'cross-corr., peak = %.3f $\mu$m' %(cc_peak/1.e4))
    
    plt.plot(lam, spec1d_B, label='Blob spectrum (467, 377)', linestyle='steps')
    
    plt.legend(loc='upper left')
    
#
def plot_F105W_backgrounds():
    """
    Check F105W backgrounds a a function of limb angle, etc.
    """
    os.chdir('/Users/brammer/WFC3/Backgrounds/BroadBand/F105W')
    filter='G141'
    alpha=0.2
    
    asns = glob.glob('*%s*orbit.png' %(filter))
    colors = np.array(['blue', 'red', 'orange'])
    #colors = np.array(['blue', 'white', 'orange'])
    fig = unicorn.plotting.plot_init(xs=6, aspect=0.7, left=0.12, right=0.03, top=0.05)
    ax = fig.add_subplot(111)
    for i, asn in enumerate(asns):
        root = asn.split('_')[0]
        print i, root
        # os.system('cat %s*%s*orbit.dat > /tmp/%s.dat' %(root[:6], filter, root))
        # bg1 = catIO.Readfile('/tmp/%s.dat' %(root), save_fits=False, force_lowercase=False)
        # bg_min = bg1.bg.min()
        # min_string = 'visit minimum'
        files=glob.glob('%s*%s*orbit.dat' %(root[:6], filter))
        bg_min = mywfc3.zodi.flt_zodi(gzfile(files[0].split('j_')[0]+'q_raw.fits'), verbose=False)
        min_string = 'zodi prediction'
        for file in files:
            bg = catIO.Readfile(file, save_fits=False, force_lowercase=False)
            cidx = bg.BrightLimb
            cidx[(bg.BrightLimb == 0) & (bg.TermAng < 120)] = 2
            #plt.scatter(bg.LimbAng, bg.bg-bg_min, c=colors[cidx], alpha=0.2, linestyle='-')
            #plt.plot(bg.LimbAng, bg.bg-bg_min, alpha=0.2, color='black')
            ### Don't or only show orbits where the limb changed during an exposure
            # limb_changed = (np.diff(bg.BrightLimb)**2).max() > 0
            # if limb_changed == 0:
            #     continue
            #
            ok = bg.BrightLimb == 1
            ax.plot(bg.LimbAng, bg.bg-bg_min, alpha=alpha/2., color='black')
            term = ~ok & (bg.TermAng > 120)
            ax.plot(bg.LimbAng[term], bg.bg[term]-bg_min, alpha=alpha, color='blue')
            ax.plot(bg.LimbAng[ok], bg.bg[ok]-bg_min, alpha=alpha, color='red')
            term = ~ok & (bg.TermAng < 120)
            ax.plot(bg.LimbAng[term], bg.bg[term]-bg_min, alpha=alpha, color='green')
    
    xli = np.arange(10, 120)

    ax.plot([-100,-100],[0,0], alpha=0.6, color='red', label='BrightLimb')
    ax.plot([-100,-100],[0,0], alpha=0.6, color='blue', label='DarkLimb')
    ax.plot([-100,-100],[0,0], alpha=0.6, color='green', label='Dark, TermAng < 120')
    #
    ax.plot(xli, 3.4564*10**(-0.06564*xli)*10, color='black', alpha=0.8, linewidth=4, label=r'STIS$\times$10 (ISR-98-21)')
    ax.plot([40,40],[-2,5], linewidth=4, linestyle='--', alpha=0.6, color='black', label='LOW-SKY; BrightLimbAng > 40')
    ax.legend(loc='upper right', prop={'size':10})
    ylims = {'F110W':(-0.5,3.5), 'F125W':(-0.2,1), 'G141':(-0.5,3.5), 'F105W':(-0.5,3.5), 'F160W':(-0.2,1)}
    ax.set_xlim(0,120); ax.set_ylim(ylims[filter])
    ax.set_xlabel('LimbAng'); ax.set_ylabel('Excess background over %s' %(min_string))
    ax.set_title(filter)
    
    unicorn.plotting.savefig(fig, '%s_backgrounds.pdf' %(filter))
    
    #########
    axr = ax.twinx()
    axr.set_ylim(np.array([-0.1,3.5])/0.65*100)
    axr.set_ylabel('%')
    
    #unicorn.plotting.savefig(fig, '../../F105W_backgrounds.pdf')
    
    #### Color bright limb by Dec.
    for asn in asns:
        print root
        root = asn.split('_')[0]
        os.system('cat %s*orbit.dat > /tmp/%s.dat' %(root[:6], root))
        bg1 = catIO.Readfile('/tmp/%s.dat' %(root), save_fits=False, force_lowercase=False)
        bg_min = bg1.bg.min()
        files=glob.glob('%s*orbit.dat' %(root[:6]))
        for file in files:
            bg = catIO.Readfile(file, save_fits=False, force_lowercase=False)
            cidx = bg.BrightLimb
            cidx[(bg.BrightLimb == 0) & (bg.TermAng < 120)] = 2
            ok = bg.BrightLimb == 1
            #plt.scatter(bg.LimbAng[ok], bg.bg[ok]-bg_min, c=np.abs(bg.DEC[ok]), vmin=0, vmax=90, alpha=0.2)
            plt.scatter(bg.read[ok], bg.LimbAng[ok], c=bg.bg[ok]-bg_min, vmin=0, vmax=4, alpha=0.2)
    #
    plt.scatter(np.arange(90), np.arange(90)*0.+4, c=np.arange(90), vmin=0, vmax=90, alpha=0.5, s=80)
    plt.xlim(0,120); plt.ylim(-1,5)
    
def geometry():
    """
    Simple plot to see how the orbit altitude compares to the surface of the earth
    """
    altitude = 550. # km
    r_earth = 6371. # km
    
    x = np.arange(-r_earth, r_earth+0.1,10)
    y = np.sqrt(r_earth**2-x**2)
    x = np.append(x,-x)
    y = np.append(y,-y)
    
    lw = 3
    plt.plot(x,y, color='green', linewidth=lw)
    plt.scatter(0,r_earth+altitude, zorder=10)
    plt.ylim(-(r_earth+5*altitude), r_earth+5*altitude)
    
    phi = np.pi/2-np.arccos(r_earth/(r_earth+altitude))
    xi, yi = np.cos(phi)*r_earth, np.sin(phi)*r_earth
    plt.plot([0,xi], [0,yi], color='blue', linewidth=lw)
    plt.plot([0,0], [0, r_earth], color='blue', linewidth=lw)
    plt.plot([0,0], [r_earth, r_earth+altitude], color='red', linewidth=lw)
    plt.plot([0, xi], [r_earth+altitude, yi], color='orange', linewidth=lw)
    
def ephem_geometry():
    import ephem
    import pyfits
    import astropy.time
    
    jit = pyfits.open('ibp321020_jit.fits.gz')
    exp = jit[2]
    spt = pyfits.getheader(gzfile(exp.header['EXPNAME'][:-1]+'q_spt.fits'), 0)
    
    spt = pyfits.open(exp.header['EXPNAME'][:-1]+'q_spt.fits.gz')
    t0 = astropy.time.Time(spt['PSTRTIME'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
    dt = astropy.time.TimeDelta(50.0, format='sec')
    
    tab = exp.data
    
    hst = ephem.Observer()
    hst.elev = 100. #550.e3
    hst.lat = tab['Latitude'][0]/180.*np.pi
    hst.lon = (tab['Longitude'][0]-360)/180*np.pi
    hst.pressure = 0.
    hst.date = ephem.Date(t0.iso)
    
    sun = ephem.Sun()
    
def compare_filters_1083():
    """
    Compare relative flux levels of the 1083 line in F098M, F105W, and F110W
    """
    import pysynphot as S
    
    hei = S.GaussianSource(1.e-16, 1.083e4, 20)
    continuum = S.FlatSpectrum(23, fluxunits='ABMag')

    for filter in ['F098M', 'F105W', 'F110W', 'F125W']:
        bp = S.ObsBandpass('wfc3,ir,%s' %(filter.lower()))
        obs_hei = S.Observation(hei, bp)
        obs_cont = S.Observation(continuum, bp)
        print '%s:  src=%0.2f  line=%.2f   linex25=%.2f  ratio=%.2f' %(filter, obs_cont.countrate(), obs_hei.countrate(), obs_hei.countrate()*25, obs_cont.countrate()/(obs_hei.countrate()*25))
        
        
    bp = S.ObsBandpass('wfc3,ir,f098m')
    obs_hei = S.Observation(hei, bp)
    obs_cont = S.Observation(continuum, bp)
    
    
    zod = S.FileSpectrum('/Users/brammer/WFC3/Backgrounds/Synphot/zodiacal_model_001.fits')
    nz = zod.renorm(SB, S.units.VegaMag, S.ObsBandpass("V"))
    bp = S.ObsBandpass('wfc3,ir,%s' %(FILTER.lower()))
    obs = S.Observation(nz, bp)
    