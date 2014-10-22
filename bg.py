"""
Get spacecraft pointing geometry from Jitter files and compare to grism backgrounds.
"""
import os

import pyfits
import numpy as np
import matplotlib.pyplot as plt

import glob

import astropy.time
import astropy.units as u

import datetime

from mywfc3.utils import gzfile
import unicorn

def cgs_to_rayleigh(flux_cgs=1.e-16, flux_cps=None, area=0.0164):
    """
    Convert CGS fluxes to Rayleighs
    
    flux has units erg/s/cm2
    area is in arcsec**2
    
    flux_cps = e-/s per WFC3 pixel in F105W
    
    # Rayleigh: 106 photons/cm2/s/Sr, or 1.58x10-11/lambda(nm) W cm2/Sr, where lambda is the wavelength of the line in nanometers.
    # 
    # http://astro.wku.edu/strolger/UNITS.txt
    
    http://www.astronomy.ohio-state.edu/~pogge/Ast871/Notes/Rayleighs.pdf    
    1 Rayleigh = 3.71546x10-14 / [wave_A] erg s-1 cm-2 arcsec-2
        
    F105W: 1.e-16 erg/s/cm2, 10A, 1.179 e/s
    """
    wavelen = 10830
    # area_Sr = area*2.35044e-11 # Area in Sr
    # flux_photon = flux_cgs/(1.986447e-8/wavelen)
    # rayleigh = 106*flux_photon/area_Sr
    
    if flux_cps is not None:
        flux_cgs = flux_cps/1.179*1.e-16
        
    rayleigh = flux_cgs/(3.715e-14/wavelen)/area
    return rayleigh
    
def go_check(im_ext='raw', force=False, asn_files=None):
    """
    Generate plots in a particular directory
    """
    import glob
    import mywfc3.bg
    
    if asn_files is None:
        asn_files = glob.glob('*jif.fits*')        
    
    for asn in asn_files:
        root = asn.split('_jif')[0]
        if (len(glob.glob('%s_*orbit.png' %(root))) == 0) | force:
            try:
                mywfc3.bg.show_orbit_limbangle(asn = [root], im_ext=im_ext)
            except ValueError:
                pass
                  
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
    asn_files = glob.glob('*jif.fits*')
    force = False
    for asn in asn_files:
        root = asn.split('_jif')[0]
        if (len(glob.glob('%s_*orbit.png' %(root))) == 0) | force:
            try:
                mywfc3.bg.show_orbit_limbangle(asn = [root])
            except ValueError:
                pass
        # try:
        #     mywfc3.bg.show_orbit_limbangle(asn = [root])
        # except:
        #     pass
            
    mywfc3.bg.show_orbit_limbangle(asn = ['ibhj44040'])
    
    ### Overall programs
    os.chdir('/Users/brammer/WFC3/GrismPrograms/CANDELS-SNe/RAW/')
    
    SKIP=True
    
    asn_files = glob.glob('*asn.fits*')

    asn_files.sort()
    for asn in asn_files:
        root = asn.split('_asn')[0]
        png = glob.glob(root+'*orbit.*')
        if (len(png) > 0) & SKIP:
            print root
            continue
        try:
            mywfc3.bg.show_orbit_limbangle(asn=[root])
        except:
            fp = open(root+'_orbit.failed','w')
            fp.write('---')
            fp.close()
            
    #### "naked" exposures
    files = glob.glob('*[a-z]j_jit.fits')
    asn = threedhst.utils.ASNFile('ic8n01020_asn.fits')
    for file in files:
        asn.exposures = [file.split('j_')[0]+'q']
        asn.product = asn.exposures[0].upper()
        root=file.split('_j')[0]
        png = glob.glob(root+'*orbit.*')
        if (len(png) > 0) & SKIP:
            print root
            continue
        #
        asn.write(root+'_asn.fits')
        mywfc3.bg.show_orbit_limbangle(asn=[root])
        
    #mywfc3.bg.show_orbit_limbangle(asn = ['ibhj20030', 'ibhj20040'])

    #### Unzip grism files in /user/brammer/WFC3_Backgrounds/GrismPrograms
    os.system('dfits *0_spt.fits |fitsort SPEC_1  |grep G1 | awk \'{print $1}\' > /tmp/log')   
    files=np.loadtxt('/tmp/log', dtype=str)
    for file in files:
        asn = threedhst.utils.ASNFile(file.replace('spt','asn'))
        for exp in asn.exposures:
            print exp
            os.system('gunzip %s*fits.gz' %(exp[:-1]))
    #
    gris_spt = np.loadtxt('/tmp/log', dtype=str)
    files=glob.glob('*asn.fits')
    asn_files = []
    for file in files:
        print file
        if file.replace('asn','spt') in gris_spt:
            asn_files.append(file)
            
      
def compute_sun_zd(root='ibffa4ckq'):
    from subprocess import Popen,PIPE
    from threedhst import catIO
    import ephem
    
    #head = pyfits.getheader(mywfc3.utils.gzfile('%s_raw.fits' %(root)))
    #h1 = pyfits.getheader(mywfc3.utils.gzfile('%s_raw.fits' %(root)), ext=1)
    #orbit = ascii.read('%sj_%s_orbit.dat' %(root[:-1], head['FILTER']))
    #ra, dec = head['RA_TARG'], head['DEC_TARG']
    
    stdout, stderr = Popen('gethead -x 0 %s_raw.fits FILTER TIME-OBS DATE-OBS' %(root), shell=True, stdout=PIPE).communicate()
    FILTER, TIME_OBS, DATE_OBS = stdout.split()
    
    orbit = catIO.Readfile('%sj_%s_orbit.dat' %(root[:-1], FILTER), save_fits=False)
    
    ### ignore time throughout exposure but account for position
    d = ephem.date('%s %s' %(DATE_OBS.replace('-','/'), TIME_OBS))
    hst = ephem.Observer()
    hst.date = d

    sun_zd = np.zeros(orbit.N)
    sun_ha = sun_zd*0.
    
    for i in range(orbit.N):
        sub_lon = orbit.longitude[i]*1
        sub_lat = orbit.latitude[i]*1
        #
        lon = sub_lon*1
        if lon > 180:
            lon = lon-360
        #
        hst.lon = '%s' %(lon)
        hst.lat = '%s' %(sub_lat)
        hst.elev = 0
        #
        sun = ephem.Sun()
        sun.compute(hst)
        sun_alt = sun.alt/np.pi*180
        sun_zd[i] = 90-sun_alt
        sun_ha[i] = (hst.sidereal_time()-sun.ra)/2/np.pi*24
        
    #print '%.3f %.3f' %(sun_alt, head['SUN_ALT'])
    return sun_zd, sun_ha
    
    if False:
        files=glob.glob('*%s*orbit.dat' %('G1'))
        for file in files:
            print file
            if os.path.exists(file.replace('dat', 'sun_zd')):
                continue
            #
            sun_zd, sun_ha = mywfc3.bg.compute_sun_zd('%sq' %(file.split('j_')[0]))
            np.savetxt(file.replace('dat', 'sun_zd'), sun_zd)
            np.savetxt(file.replace('dat', 'sun_ha'), sun_ha)
            
def compute_hst_azimuth(sub_lat, sub_lon, root='ibffa4ckq'):
    """
    Compute azimuth for an HST observation, assuming the orbital 
    parameters have already been extracted from the jitter files.
    """
    import mywfc3.utils
    from astropy.io import ascii
    import ephem
    
    head = pyfits.getheader(mywfc3.utils.gzfile('%s_raw.fits' %(root)))
    #orbit = ascii.read('%sj_%s_orbit.dat' %(root[:-1], head['FILTER']))
    ra, dec = head['RA_TARG'], head['DEC_TARG']
    
    d = ephem.date('%s %s' %(head['DATE-OBS'].replace('-','/'), head['TIME-OBS']))
    hst = ephem.Observer()
    hst.date = d
    #lon = orbit['Longitude'][0]
    lon = sub_lon*1
    if lon > 180:
        lon = lon-360
    #
    hst.lon = '%s' %(lon)
    hst.lat = '%s' %(sub_lat)
    hst.elev = 600.e3
        
    targ = ephem.star('Arcturus')
    targ._ra, targ._dec = ra/180*np.pi, dec/180.*np.pi
    targ._pmdec, targ._pmra, = 0, 0
    targ.compute(hst)
    
    sun = ephem.Sun()
    sun.compute(hst)
    sun_zd = 90-sun.alt/np.pi*180
    
    az = targ.az/np.pi*180
    return az
    
def compute_shadow(root='ib5x15vlq', save_fits=False):
    
    h = pyfits.getheader(root+'_raw.fits')
    f = h['FILTER']
    orbit = catIO.Readfile('%sj_%s_orbit.dat' %(root[:-1], f), save_fits=False)
    asn = h['ASN_ID'].strip().lower()
    jif = pyfits.open(gzfile('%s_jif.fits' %(asn)))
    for ext in jif[1:]:
        if ext.header['EXPNAME'] == root[:-1]+'j':
            pstr = astropy.time.Time(ext.header['STARTIME'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
            shadow_in = astropy.time.Time(ext.header['SHADOENT'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
            shadow_out = astropy.time.Time(ext.header['SHADOEXT'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
            dt_in = (shadow_in-pstr).sec
            dt_out = (shadow_out-pstr).sec
            in_shadow = ((orbit.seconds <= dt_out) & (dt_out <= orbit.seconds.max())) | ((orbit.seconds >= dt_in) & (dt_in > 0))
            #orbit = catIO.Readfile('%sj_%s_orbit.dat' %(root[:-1], f), save_fits=False)
            orbit.addColumn('in_shadow', in_shadow*1)
            orbit.write_text(orbit.filename)
            
def show_orbit_limbangle(asn = ['ib3701050', 'ib3701060'], ymax=3.8, im_ext='raw', tstr=None):
    
    import scipy.ndimage as nd
    import mywfc3.utils
    import mywfc3.zodi
    import mywfc3.etc_zodi
    
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
        spt_file = glob.glob('%s*spt.fits*' %(ext.header['EXPNAME'][:-1]))[0]
        expname = spt_file.split('_spt')[0]
        #expname = ext.header['EXPNAME'][:-1]+'q'
        print expname
        spt = pyfits.getheader(spt_file, 0)
        if i == 0:
            targname = spt['TARGNAME']
            ax3.text(0.5, 0.95, targname, ha='center', va='top', transform=ax3.transAxes)
            try:
                zodi_obj = mywfc3.etc_zodi.ZodiacalLight(mywfc3.utils.gzfile(expname+'_raw.fits'))
                zodi_obj.eval_filter(verbose=True)
                zodi = zodi_obj.countrate
            except:
                zodi = mywfc3.zodi.flt_zodi(mywfc3.utils.gzfile(expname+'_raw.fits'), verbose=False)                
        #
        #### Start/stop times
        pstr = astropy.time.Time(spt['PSTRTIME'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
        pstp = astropy.time.Time(spt['PSTPTIME'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
        if i == 0:
           if tstr is None:
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
        #print gzfile('%s_%s.fits' %(expname, im_ext)), '%s_%s.fits' %(expname, im_ext)
        ima = pyfits.open(gzfile('%s_%s.fits' %(expname, im_ext)))
        if im_ext == 'ima':
            GAIN = 1
        else:
            GAIN = 2.5
        #
        FILTER = ima[0].header['FILTER']
        time, ramp, reads = get_bg_ramp(ima)
        dt = np.diff(time)
        ok = dt > 24
        ymax = np.maximum(ymax, (ramp/dt*GAIN).max()*1.05)
        ax3.plot(((pstr-tstr).sec + time[1:][ok])/60., (ramp/dt*GAIN)[ok], color=colors[i], linewidth=2)
        ax3.plot(((pstr-tstr).sec + time[1:][ok])/60., ramp[ok]*0+zodi, color='black', linestyle='--', alpha=0.5)
        #
        ##  Get shadow information from the JIF file
        ## JIF Shadow
        shadow_flag = ext.data['Seconds']*0
        if jif is not None:
            shadow_in = astropy.time.Time(jif[i].header['SHADOENT'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
            shadow_out = astropy.time.Time(jif[i].header['SHADOEXT'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
            dt_in = (shadow_in-pstr).sec
            dt_out = (shadow_out-pstr).sec
            print 'JIF: in=%.3f out=%.3f' %(dt_in, dt_out)
            if dt_in > 0:
                yi = np.interp(dt_in, (time[1:][ok]), (ramp/dt*GAIN)[ok])
                ax3.scatter(((pstr-tstr).sec + dt_in)/60., yi, marker='o', s=40, color='green', alpha=0.8)
            #
            if (dt_out > 0) & (dt_out < time.max()):
                yi = np.interp(dt_out, (time[1:][ok]), (ramp/dt*GAIN)[ok])
                ax3.scatter(((pstr-tstr).sec + dt_out)/60., yi, marker='o', s=40, color='red', alpha=0.8)
            #
            in_shadow = ext.data['Seconds'] < -1
            if (dt_in < 0) & (dt_out < 0):
                if dt_out < dt_in:
                    in_shadow = (ext.data['Seconds'] > -1)
            #
            if (dt_in > 0):
                in_shadow = (ext.data['Seconds'] >= dt_in)
            #
            if (dt_out > 0):
                in_shadow = (ext.data['Seconds'] <= dt_out)
            #
            #print in_shadow.sum(), dt_in, dt_out
            #ax3.plot(((pstr-tstr).sec + ext.data['Seconds'][in_shadow])/60., ext.data['Seconds'][in_shadow]*0.+2, color='green')
            shadow_flag[in_shadow] = 1
            shadow_flag[~in_shadow] = 0
            
        #### Make log
        idx = np.arange(len(time)-1)[ok]
        fp = open('%s_%s_orbit.dat' %(ext.header['EXPNAME'], FILTER), 'w')
        dtypes = {'E':'%.3f', 'D':'%.3f', 'I':'%d'}
        line = '# name read  bg minbg zodi shadow'
        for c in ext.data.columns:
            line += ' %s' %(c.name)
        fp.write(line+'\n')
        
        for il in idx:
            line = '%s_%02d %2d %.3f %.3f %.3f' %(ext.header['EXPNAME'], il, il, ramp[il]/dt[il]*GAIN, np.min((ramp/dt*GAIN)[ok]), zodi)
            #
            shadow = np.interp(time[1:][il]-dt[il]/2., ext.data['Seconds'], shadow_flag)
            line += ' %d' %(shadow)
            for c in ext.data.columns:
                val = np.interp(time[1:][il]-dt[il]/2., ext.data['Seconds'], ext.data[c.name])
                line += ' ' + dtypes[c.format] %(val) # %(str(np.cast[dtypes[c.format]](val)))
            #
            fp.write(line+'\n')
        fp.close()
        
        #img = nd.convolve(flt[1].data, np.ones((2,2))/4.)/med_background
        last_read = (ima['SCI', 1].data-ima['SCI',len(time)].data)[5:-5, 5:-5]*GAIN
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
    ax1.set_ylabel('LimbAng')
    ax1.set_ylim(0,99)
    xl = ax1.get_xlim()
    ax1.plot([-10,xl[1]], [20,20], linestyle=':')
    ax1.set_xlim(-10, xl[1])
    ax3.set_xlim(-10, xl[1])
    #ax3.legend(loc='upper left', prop={'size':8})
    
    ax3.set_ylim(0,ymax)
    ax3.set_xticklabels([])
    ax3.set_ylabel('backg (electrons/s)')
    ax2 = fig.add_axes((0.66,0.14,0.33,0.4))
    
    
    map = init_map(ax=ax2)
    for i, ext in enumerate(jit):
        draw_map_latlon(map, ext.data['Latitude'], ext.data['Longitude'], color=colors[i], linewidth=2, ext=ext)
        
    tt = np.cast[int](tstr.iso.replace('-',' ').replace(':',' ').replace('.', ' ').split())
    
    #date = datetime.utcnow()
    date = datetime.datetime(tt[0], tt[1], tt[2], tt[3], tt[4], tt[5], tt[6])
    CS=map.nightshade(date, alpha=0.2, color='black')
    ax2.set_title(tstr.iso)
    
    #plt.savefig('%s_%s_orbit.png' %(asn[0], FILTER))
    unicorn.plotting.savefig(fig, '%s_%s_orbit.png' %(asn[0], FILTER))
    #plt.close()
    
    return tstr
    
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
    
def draw_map_latlon(map, lat, lon, ext=None, *args, **kwargs):
    w = lon < 180
    e = lon >= 180
    xpt, ypt = map(lon[w], lat[w])
    map.plot(xpt, ypt, alpha=0.5, **kwargs)
    xpt, ypt = map(lon[e], lat[e])
    map.plot(xpt, ypt, alpha=0.5, **kwargs)
    xpt, ypt = map(lon[0], lat[0])
    map.scatter(xpt, ypt, alpha=0.5, **kwargs)
    if ext is not None:
        root = ext.header['EXPNAME'][:-1]+'q'
        az = compute_hst_azimuth(lat[0], lon[0], root=root)
        #print az
        dx = np.cos((90-az)/180*np.pi)
        dy = np.sin((90-az)/180*np.pi)
        xpt2, ypt2 = map(lon[0]+dx*10, lat[0]+dy*10)
        map.plot([xpt, xpt2], [ypt, ypt2], color='red')
        
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
        ds9.view((sky[i,:,:]))
        
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
    
def split_multiaccum(ima, scale_flat=True):
    
    if scale_flat:
        FLAT_F140W = pyfits.open(os.path.join(os.getenv('iref'), 'uc721143i_pfl.fits'))[1].data
    else:
        FLAT_F140W = 1.
    
    if ('ima' in ima.filename()) & (ima[0].header['FLATCORR'] == 'COMPLETE'):
        FLAT_F140W = 1
        
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

def evaluate_bg_sn(raw='ib3712lyq_raw.fits', npix=4):
    """
    Plot S/N as a function of time for different source counts
    """
    from threedhst import catIO
    import glob
    import numpy as np
    import matplotlib.pyplot as plt
    
    s = 1 # e/s
    bg = catIO.Readfile(glob.glob('%sj_*_orbit.dat' %(raw.split('q_')[0]))[0], save_fits=False)
    
    dt = np.ones(bg.N)*100
    t = np.cumsum(dt)
    RN = 20
    
    signal = np.cumsum(s*dt)
    var = npix*RN**2 + np.cumsum((s+npix*bg.bg)*dt)
    var_zodi = npix*RN**2 + np.cumsum((s+npix*bg.zodi)*dt)
    SN = signal/np.sqrt(var)
    SN_zodi = signal/np.sqrt(var_zodi)
    plt.plot(t, SN/SN[3], label=raw.split('_')[0])
    plt.plot(t, SN_zodi/SN[3], color='black', linestyle='--')
        
    print np.ptp(bg.bg[1:])
    
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
        #bg_min = mywfc3.zodi.flt_zodi(gzfile(files[0].split('j_')[0]+'q_raw.fits'), verbose=False)
        min_string = 'zodi prediction'
        for file in files:
            bg = catIO.Readfile(file, save_fits=False, force_lowercase=False)
            bg_min = bg.zodi
            #
            ok = bg.shadow == 1
            ax.plot(bg.LimbAng, bg.bg-bg_min, alpha=alpha/2., color='black')
            ax.plot(bg.LimbAng[ok], (bg.bg-bg_min)[ok], alpha=alpha, color='blue')
            ax.plot(bg.LimbAng[~ok & (bg.BrightLimb == 1)], (bg.bg-bg_min)[~ok & (bg.BrightLimb == 1)], alpha=alpha, color='red')
            ax.plot(bg.LimbAng[~ok & ~(bg.BrightLimb == 1)], (bg.bg-bg_min)[~ok & ~(bg.BrightLimb == 1)], alpha=alpha, color='green')
    
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
    
    #
    bg = catIO.Readfile('master.dat', force_lowercase=False, save_fits=False)
    s = bg.shadow == 1
    b = bg.BrightLimb == 1
    
    plt.scatter(bg.LimbAng[~s], (bg.bg-bg.zodi)[~s], alpha=0.1)
    
    sun = np.zeros(bg.N)
    for i in range(bg.N):
        sun[i] = sun_angle['%sq' %(bg.name[i].split('j_')[0])]
    #
    vmax=140
    
    plt.scatter(bg.LimbAng[~s], (bg.bg/bg.zodi)[~s], alpha=0.5, c=sun[~s], vmin=45, vmax=vmax)
    plt.ylim(0,10); plt.xlim(0,120); plt.ylabel('BG / zodi')
    cb = plt.colorbar(); cb.set_label('Sun angle')
    plt.savefig('ShadowFalse.png')
    
    plt.scatter(bg.LimbAng[s], (bg.bg/bg.zodi)[s], alpha=0.5, c=sun[s], vmin=45, vmax=vmax)
    plt.ylim(0,10); plt.xlim(0,120); plt.ylabel('BG / zodi')
    cb = plt.colorbar(); cb.set_label('Sun angle')
    plt.savefig('ShadowTrue.png')

def show_with_f10(log=False):
    """ 
    Make plots of LimbAngle vs background excess colored by Solar 10.7cm Flux
    """   
    import mywfc3.bg
    from threedhst import catIO
    import cPickle as pickle
    import scipy
    
    ###### Figure of tracks colored by sun angle
    cmap = 'cool'
    cmap = 'RdYlGn_r'
    
    la, ra, ta, ba = 0.06, 0.1, 0.03, 0.13
    fig = unicorn.plotting.plot_init(xs=8, aspect=0.4, left=0.01, right=0.01, top=0.01, bottom=0.01, wspace=0, hspace=0, NO_GUI=True, use_tex=True)
        
    axes=[]
    dx = (1-la-ra)/3.
    for i in range(3):
        ax = fig.add_axes((la+dx*i, ba, dx, 1-ba-ta))
        axes.append(ax)
    #
    cax = fig.add_axes(((1-ra)+ra/10., ba, ra/4., 1-ta-ba))
    
    filters = ['F105W','G102','G141']
    #filters = ['F105W']
    os.chdir('/user/brammer/WFC3_Backgrounds/F105W')
    #filters = ['F125W','G102','G141']
    #os.chdir('/user/brammer/WFC3_Backgrounds/F125W')
    for i in range(len(filters)):
        if filters[i].startswith('G'):
            os.chdir('../GrismPrograms')
        #
        ax = fig.axes[i]
        mywfc3.bg.single_sunangle(filter=filters[i], alpha=0.2, log=log, ax=ax, cax=cax, cmap=cmap)
        ax.text(0.95, 0.95, filters[i], size=13, ha='right', va='top', transform=ax.transAxes)
        
    #
    Clabel, vc = 'Solar 10.7 cm flux (sfu)', (70,160)
    s = ax.scatter(np.ones(2),np.ones(2)*1000, c=vc, vmin=vc[0], vmax=vc[1], alpha=0.5, cmap=cmap)
    cb = fig.colorbar(s, cax=cax)
    cb.set_label(Clabel)
    
    for i in range(1,3):
        fig.axes[i].set_ylabel('')
        fig.axes[i].set_yticklabels([])
        
    unicorn.plotting.savefig(fig, '../excess_solar_F10x.pdf', increment=False)
    
def sunangle(filter='F105W', alpha=0.2):
    import mywfc3.bg
    from threedhst import catIO
    import cPickle as pickle
    import scipy
    
    asns = glob.glob('*%s*orbit.png' %(filter))
    fp = open('f105w_sun_coords.pkl','rb')    
    sun_angle = pickle.load(fp)
    fp.close()
    
    ###### Figure of tracks colored by sun angle
    la, ra, ta, ba = 0.06, 0.1, 0.03, 0.13
    fig = unicorn.plotting.plot_init(xs=8, aspect=0.4, left=0.01, right=0.01, top=0.01, bottom=0.01, wspace=0, hspace=0, NO_GUI=True, use_tex=True)
    
    master = {131:{'label':'Shadow', 's':1, 'b':0},
              132:{'label':'In sunlight, dark limb', 's':0, 'b':0},
              133:{'label':'In sunlight, bright limb', 's':0, 'b':1}}
              
    axes = []
    dx = (1-la-ra)/3.
    for i, key in enumerate(master.keys()):
        ax = fig.add_axes((la+dx*i, ba, dx, 1-ba-ta))
        axes.append(ax)
    
    for i, asn in enumerate(asns):
        root = asn.split('_')[0]
        print i, root
        files=glob.glob('%s*%s*orbit.dat' %(root[:6], filter))
        min_string = 'zodi prediction'
        for file in files:
            bg = catIO.Readfile(file, save_fits=False, force_lowercase=False)
            bg_min = bg.zodi
            ### Correction for shadow zodi as a function of sun angle
            sun_ang = sun_angle['%sq' %(bg.name[0].split('j_')[0])]
            p = np.array([ -3.69229110e-10,   2.39078911e-07,  -5.86843639e-05, 6.75812691e-03,  -3.62093498e-01,   8.16467389e+00])
            bg_min *= scipy.polyval(p, sun_ang)
            #
            ok = bg.shadow == 1
            for j, key in enumerate(master.keys()):
                ax = axes[j]
                idx = (bg.shadow == master[key]['s']) & (bg.BrightLimb == master[key]['b'])
                yFlux, Ylabel = (bg.bg/bg_min), 'Background / zodi prediction'
                #yFlux, Ylabel = (bg.bg - bg_min), 'Background - zodi prediction'
                sun_ha = np.loadtxt(file.replace('dat','sun_ha')) #, 'SUN ZD'
                #
                xAngle, Alabel, xr = bg.LimbAng, 'LimbAng', (10,119)
                #xAngle, Alabel, xr = bg.LOS_Zenith, 'LOS Zenith', (1,100)
                #xAngle, Alabel, xr = np.loadtxt(file.replace('dat','sun_zd')), 'SUN ZD', (1,178)
                xAngle, Alabel, xr = sun_ha, 'Solar HA', (-14, 14)
                #idx = idx & (sun_ha < -2)
                ax.plot(xAngle, yFlux, color='black', alpha=0.03)
                if idx.sum() > 0:
                    #color_index, Clabel, vc = np.ones(bg.N)[idx]*sun_ang, 'Target-Sun Angle', (45, 140)
                    color_index, Clabel, vc = bg.LOS_Zenith[idx], 'LOS Zenith', (1,100)
                    ## color_index, Clabel, vc = np.loadtxt(file.replace('dat','sun_zd'))[idx], 'SUN ZD', (45,140)
                    #color_index, Clabel, vc = sun_ha[idx], 'Solar Hour Angle', (-5,5)
                    
                    mywfc3.bg.color_lines(xAngle[idx], yFlux[idx], color_index, alpha=alpha, linewidth=1.2, ax=ax, vmin=vc[0], vmax=vc[1], cmap='jet', clip_dx=5)
    
    for j, key in enumerate(master.keys()):
        axes[j].plot([0,200], np.ones(2)*('/' in Ylabel), color='white', linewidth=3, alpha=0.8)
        axes[j].plot([0,200], np.ones(2)*('/' in Ylabel), color='black', linewidth=1.3, alpha=0.8, linestyle='--')
        axes[j].set_xlim(xr); axes[j].set_ylim(-0.5,5)
        axes[j].set_xlabel(Alabel)
        axes[j].text(0.5, 0.95, master[key]['label'], ha='center', va='top', transform=axes[j].transAxes)
            
    #
    s = axes[0].scatter(np.ones(2),np.ones(2)*1000, c=vc, vmin=vc[0], vmax=vc[1], alpha=0.5)
    cax = fig.add_axes(((1-ra)+ra/10., ba, ra/4., 1-ta-ba))
    cb = fig.colorbar(s, cax=cax)
    cb.set_label(Clabel)
    
    ### Ylog
    if '/' in Ylabel:
        for j, key in enumerate(master.keys()):
            axes[j].semilogy()
            axes[j].set_ylim(0.7,15)
        
        axes[0].set_yticklabels([1,10])
        axes[0].set_yticks([1,10])
        
    axes[0].set_ylabel(Ylabel)
    for i in range(1,3):
        axes[i].set_yticklabels([])
    
    unicorn.plotting.savefig(fig, 'F105W_background_SunAngle.pdf', increment=True)
    
    if False:
        #### Look into offset of shadow zodi / obs as a function of sun angle
        bgf = catIO.Readfile('master.dat', force_lowercase=False, save_fits=False)
        bgf.exposure = bgf.name
        for i in range(bgf.N):
            bgf.exposure[i] = bgf.name[i].split('j_')[0]+'q'
        #
        # solar flux
        sun_f10 = np.loadtxt('master.sun_f10')
        
        sun_zd = np.loadtxt('master.sun_zd')
        sun_ha = np.loadtxt('master.sun_ha')
        sun_ha[sun_ha < -12] += 24 #-sun_ha[sun_ha < -12]-24
        sun_ha[sun_ha > 12] -= 24 #-sun_ha[sun_ha > -12]+24
        sun_a = np.zeros(bgf.N)
        for i in range(bgf.N):
            print i
            sun_a[i] = sun_angle['%sq' %(bgf.name[i].split('j_')[0])]
        
        #
        shadow = sun_zd > 130
        plt.scatter(sun_a[shadow], (bgf.bg/bgf.zodi)[shadow], alpha=0.1)
        p = scipy.polyfit(sun_a[shadow], (bgf.bg/bgf.zodi)[shadow], 5)
        xarr = np.arange(0,180)
        plt.plot(xarr, scipy.polyval(p, xarr))
        #
        ### Check for difference with morning and twilight (hour angle)
        in_sun = sun_zd <  100
        plt.scatter(sun_zd[in_sun], (bgf.bg/bgf.zodi)[in_sun], c=sun_a[in_sun], vmin=45, vmax=140, alpha=0.1)

        plt.scatter(sun_zd[in_sun], (bgf.bg/bgf.zodi)[in_sun], c=sun_ha[in_sun], vmin=-5, vmax=5, alpha=0.4)
        plt.colorbar()

        ####
        from subprocess import Popen,PIPE
        flux10 = catIO.Readfile('../Solar_10.7_fluxtable.txt')
        bgf.f10 = bgf.bg*0
        files = glob.glob('i*orbit.dat')
        for file in files:
            print file
            exp = file.split('j_')[0]
            stdout, stderr = Popen('gethead -x 0 %sq_raw.fits EXPSTART NSAMP' %(exp), shell=True, stdout=PIPE).communicate()
            mjd, nsamp = np.cast[float](stdout.split())
            jd0 = mjd+24.e5+0.5
            f10_exp = np.interp(jd0, flux10.fluxjulian, flux10.fluxobsflux)
            np.savetxt(file.replace('orbit.dat', 'orbit.sun_f10'), np.ones(nsamp-1).T*f10_exp)
            
def single_sunangle(filter='F105W', alpha=0.2, log=False, ax=None, cax=None, cmap='jet'):
    import mywfc3.bg
    from threedhst import catIO
    import cPickle as pickle
    import scipy
    
    asns = glob.glob('*%s*orbit.png' %(filter))
    asns.sort()
    fp = open('f105w_sun_coords.pkl','rb')    
    sun_angle = pickle.load(fp)
    fp.close()
    
    ###### Figure of tracks colored by sun angle
    aspect, la, ra, ta, ba = 0.7, 0.06, 0.11, 0.02, 0.08
    aspect, la, ra, ta, ba = 0.5, 0.06, 0.105, 0.02, 0.105

    fig = None
    if ax is None:
        fig = unicorn.plotting.plot_init(xs=7, aspect=aspect, left=0.01, right=0.01, top=0.01, bottom=0.01, wspace=0, hspace=0, NO_GUI=True, use_tex=True)
        ax = fig.add_axes((la, ba, 1-la-ra, 1-ba-ta))
    
    p = {}
    
    p['F105W'] = np.array([ -3.69229110e-10,   2.39078911e-07,  -5.86843639e-05, 6.75812691e-03,  -3.62093498e-01,   8.16467389e+00])
    p['F125W'] = np.array([ -3.69229110e-10,   2.39078911e-07,  -5.86843639e-05, 6.75812691e-03,  -3.62093498e-01,   8.16467389e+00])
    
    #### M31 WISP
    bad_visits = ['ibtt69', 'ibtt70', 'ibtt71', 'ibtt72', 'ibtt73', 'ibtt74', 'ibtt75', 'ibtt76', 'ibtt77', 'ibtt78', 'ibtt79', 'ibtt80', 'ibtt81', 'ibtt82', 'ibtt83', 'ibtt84', 'ibtt98', 'ibtt1s', 'ibtt1t', 'ibtt1u', 'ib8ceo', 'ib8cep']
    
    #### WISP extremely bright star (H=7)
    bad_visits.extend(['ic3t85', 'ic3t86', 'ic3t87', 'ic3t88', 'ic3t89'])
    for i, asn in enumerate(asns[:]):
        root = asn.split('_')[0]
        if root[:6] in bad_visits:
            continue
        #
        print i, root
        files=glob.glob('%s*%s*orbit.dat' %(root[:6], filter))
        min_string = 'zodi prediction'
        for file in files:
            bg = catIO.Readfile(file, save_fits=False, force_lowercase=False)
            bg_min = bg.zodi
            ### Correction for shadow zodi as a function of sun angle
            sun_ang = sun_angle['%sq' %(bg.name[0].split('j_')[0])]
            #bg_min *= scipy.polyval(p[filter], sun_ang)
            #
            yFlux, Ylabel = (bg.bg - bg_min), 'Background - zodi prediction, e-/s'
            if log:
                yFlux, Ylabel = (bg.bg/bg_min), 'Background / zodi prediction'
        
            sun_ha = np.loadtxt(file.replace('dat','sun_ha')) #, 'SUN ZD'
            sun_ha[sun_ha < -12] += 24 #-sun_ha[sun_ha < -12]-24
            sun_ha[sun_ha > 12] -= 24 #-sun_ha[sun_ha > -12]+24
            #
            # print file, np.diff(sun_ha).min(), np.diff(sun_ha).max()   
            #xAngle, Alabel, xr = bg.LimbAng, 'LimbAng', (10,119)
            #ok = bg.shadow == 0
            #xAngle, Alabel, xr = bg.LOS_Zenith, 'LOS Zenith', (1,100)
            #xAngle, Alabel, xr = np.loadtxt(file.replace('dat','sun_zd')), 'SUN ZD', (0,178)
            xAngle, Alabel, xr = sun_ha, 'Subpoint Solar HA (hours)', (-12,12)
            ok = xAngle > -1000
            #### Remove some dense clusters
            if filter == 'F125W':
               ok = ok & (yFlux < 1.6*(1+0.56*('d / z' in Ylabel)))
            #ok = (bg.LimbAng >= 40) & (bg.shadow == 1)
            #
            #ax.plot(xAngle, yFlux, color='black', alpha=0.03)
            #print ok.sum()
            
            #
            ### Colored line segment
            #color_index, Clabel, vc = np.ones(bg.N)*sun_ang, 'Target-Sun Angle', (45, 140)
            #color_index, Clabel, vc = bg.LOS_Zenith, 'LOS Zenith', (1,100)
            color_index, Clabel, vc = bg.LimbAng, 'Limb Angle', (15,100)
            ## color_index, Clabel, vc = np.loadtxt(file.replace('dat','sun_zd')), 'SUN ZD', (45,140)
            #color_index, Clabel, vc = sun_ha, 'Solar Hour Angle', (-5,5)
            sun_f10 = np.loadtxt(file.replace('dat','sun_f10')) #, 'SUN ZD'
            color_index, Clabel, vc = sun_f10, 'Solar 10.7cm flux', (70,160)
            xAngle, Alabel, xr = bg.LimbAng, 'LimbAng', (10,119)
            
            ### PMCC's "airmass" idea
            xAngle, Alabel, xr = 1/np.cos((104-bg.LimbAng)/180*np.pi), 'LimbAM', (1,6.9)
            
            ok = ok & (np.abs(sun_ha) < 5)
            
            if ok.sum() == 0:
                continue
                        
            mywfc3.bg.color_lines(xAngle[ok], yFlux[ok], yFlux[ok]*0, alpha=0.05, linewidth=1.2, ax=ax, vmin=0, vmax=1, cmap='gray', clip_dx=5+100*('hours' not in Alabel))

            mywfc3.bg.color_lines(xAngle[ok], yFlux[ok], color_index[ok], alpha=alpha, linewidth=1.2, ax=ax, vmin=vc[0], vmax=vc[1], cmap=cmap, clip_dx=5+100*('hours' not in Alabel))
    
    ax.plot([-200,200], np.ones(2)*('d / z' in Ylabel), color='white', linewidth=3, alpha=0.8)
    ax.plot([-200,200], np.ones(2)*('d / z' in Ylabel), color='black', linewidth=1.3, alpha=0.8, linestyle='--')
    ax.set_xlim(xr); ax.set_ylim(-0.5,6)
    ax.set_xlabel(Alabel)
                
    if fig is not None:
        cax = fig.add_axes(((1-ra)+ra/10., ba, ra/4., 1-ta-ba))
        s = ax.scatter(np.ones(2),np.ones(2)*1000, c=vc, vmin=vc[0], vmax=vc[1], alpha=0.5)
        cb = fig.colorbar(s, cax=cax)
        cb.set_label(Clabel)
    
    yshade = [5,20,5.2,5.6]
    ### Ylog
    if 'd / z' in Ylabel:
        ax.semilogy()
        ax.set_ylim(0.7,15)
        
        ax.set_yticklabels([1,2,5,10])
        ax.set_yticks([1,2,5,10])
        yshade = [9,20,10,12]
        
    ax.set_ylabel(Ylabel) # + ', %s' %(filter))
    if fig is not None:
        ax.text(0.1, 0.7, filter, size=16, ha='left', va='bottom', transform=ax.transAxes)
    
    if 'HA' in Alabel:
        ax.fill_between([-14,-8], np.ones(2)*yshade[0], np.ones(2)*yshade[1], color='0.2', alpha=0.5)
        ax.text(-9, yshade[2], 'Shadow', ha='right', va='bottom')
        ax.fill_between([8,14], np.ones(2)*yshade[0], np.ones(2)*yshade[1], color='0.2', alpha=0.5)
        ax.text(9, yshade[2], 'Shadow', ha='left', va='bottom')
        
        ax.fill_between([-8,-5], np.ones(2)*yshade[0], np.ones(2)*yshade[1], color='0.6', alpha=0.5)
        ax.text(-5.5, yshade[2], 'Dawn', ha='right', va='bottom')
        ax.fill_between([5,8], np.ones(2)*yshade[0], np.ones(2)*yshade[1], color='0.6', alpha=0.5)
        ax.text(5.5, yshade[2], 'Dusk', ha='left', va='bottom')

        ax.text(0, yshade[2], 'Day', ha='center', va='bottom')
        ax.text(0, yshade[3], r'\textit{HST} in ', ha='center', va='bottom')
        
    if fig is not None:
        unicorn.plotting.savefig(fig, filter+'_background_SunAngle_single.pdf', increment=True)
    
       
def color_lines(x, y, z, vmin=0, vmax=1, alpha=1, linewidth=2, cmap='spectral', clip_dx=None, ax=None):
    """
    Color line segments x vs y using the array `z`, normalized to (vmin,vmax)
    
    http://wiki.scipy.org/Cookbook/Matplotlib/MulticoloredLine
    """
    import matplotlib.collections
    
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    if clip_dx is not None:
        dx = segments[:,1,0]-segments[:,0,0]
        clip = np.abs(dx) < clip_dx
        if clip.sum() == 0:
            return False
        
        segments = segments[clip,:,:]
        
    lc = matplotlib.collections.LineCollection(segments, cmap=plt.get_cmap(cmap), norm=plt.Normalize(vmin, vmax))
    lc.set_array(z)
    lc.set_linewidth(linewidth)
    lc.set_alpha(alpha)
    if ax is None:
        ax = plt.gca()
    #
    ax.add_collection(lc)
    
#### Sun angle
def get_sun_angle():
    """
    Get SUN coords from SPT file
    
    dfits *q_spt.fits |fitsort PSTRTIME UTC0 RA_TARG DEC_TARG RA_SUN DEC_SUN > sun_coords.dat
    """
    import pyfits
    import astropy.coordinates as c
    import astropy.units as u
    import cPickle as pickle
    
    if hasattr(c, 'ICRS'):
        icrs = c.ICRS
    else:
        icrs = c.ICRSCoordinates
        
    #h = pyfits.getheader(spt)
    
    sun = catIO.Readfile('sun_coords.dat', save_fits=False, force_lowercase=False)
    sun_angle = {}
    
    for i in range(sun.N):
        print i
        key = sun.FILE[i].split('_spt')[0]
        ra, dec = sun.RA_TARG[i], sun.DEC_TARG[i]
        ra_sun, dec_sun = sun.RA_SUN[i], sun.DEC_SUN[i]
        c_targ = icrs(ra, dec, unit=(u.deg, u.deg))
        c_sun = icrs(ra_sun, dec_sun, unit=(u.deg, u.deg))
        sun_angle[key] = c_targ.separation(c_sun).degrees #value
    #
    fp = open('f105w_sun_coords.pkl','wb')
    pickle.dump(sun_angle, fp)
    fp.close()

    fp = open('f105w_sun_coords.pkl','rb')    
    sun_angle = pickle.load(fp)
    fp.close()
    
    
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
    
    #### Height of shadow at horizon and zenith as a function of solar depression angle
    #### from the earth surface
    alpha = np.arange(0,90.,1.) # 90 - depression angle
    theta = (180-(90-alpha))/2.
    h_horizon = r_earth * (1./np.sin(theta/180.*np.pi)-1)
    h_zenith = r_earth * (1./np.sin(alpha/180.*np.pi)-1)
    
    plt.plot(90-alpha, h_horizon, label='h_horiz')
    plt.plot(90-alpha, h_zenith, label='h_zenith')
    plt.plot([0,90], np.ones(2)*altitude, label='HST')
    plt.semilogy()
    plt.legend(loc='upper left')
    
    
    
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
    
    line_level = 1 # e / s / pix
    hei = S.GaussianSource(1.67e-16/2.0*line_level, 1.083e4, 20)
    #hei = S.GaussianSource(8.403e-17, 1.083e4, 20)
    continuum = S.FlatSpectrum(23, fluxunits='ABMag')
    #continuum = S.FlatSpectrum(23*1.e-9, fluxunits='photlam')

    zodi = S.FileSpectrum('/Users/brammer/WFC3/Backgrounds/Synphot/zodiacal_model_001.fits')
    nz = zodi.renorm(25, S.units.VegaMag, S.ObsBandpass("V"))
    bp = S.ObsBandpass('wfc3,ir,%s' %('F105W'.lower()))
    obs_zodi = S.Observation(nz, bp)
    zodi_refmag = 25+2.5*np.log10(2.84/30.63)
    zodi_level = 1 # e/s/pix
    refmag = 25+2.5*np.log10(obs_zodi.countrate()/zodi_level)
    nz = zodi.renorm(refmag, S.units.VegaMag, S.ObsBandpass("V"))
    
    #wfc3_ir_dark = 0.0192
    thermal = 1.556/30.63
    dark = 1.532/30.63
    RN = 14.6
    
    for filter in ['F098M', 'F105W', 'F110W', 'F125W', 'F140W', 'F160W', 'G102','G141']:
        bp = S.ObsBandpass('wfc3,ir,%s' %(filter.lower()))
        obs_hei = S.Observation(hei, bp)
        obs_cont = S.Observation(continuum, bp)
        obs_zodi = S.Observation(nz, bp)
        print '%s:  src=%0.2f  line=%.2f   linex25=%.2f  ratio=%.2f zodi=%.2f' %(filter, obs_cont.countrate(), obs_hei.countrate(), obs_hei.countrate()*25, obs_cont.countrate()/(obs_hei.countrate()*25), obs_zodi.countrate())

    #### Check optical filters -> nothing
    # for bp_str in ['acs,wfc1,f814w', 'acs,wfc1,f850lp', 'wfc3,uvis1,f850lp']:
    #     bp = S.ObsBandpass(bp_str)
    #     obs_hei = S.Observation(hei, bp)
    #     obs_cont = S.Observation(continuum, bp)
    #     obs_zodi = S.Observation(nz, bp)
    #     print '%s:  src=%0.2f  line=%.2f   linex25=%.2f  ratio=%.2f zodi=%.2f' %(bp, obs_cont.countrate(), obs_hei.countrate(), obs_hei.countrate()*25, obs_cont.countrate()/(obs_hei.countrate()*25), obs_zodi.countrate())
       
    
    line_flam = 1.e-16/1.19
    pix = 0.128254
    flux_rayleigh = line_flam/pix**2/(3.7e-14/1.083e4)
    
    ### G102
    line_flam = 1.e-16/1.00
    pix = 0.128254
    flux_rayleigh = line_flam/pix**2/(3.7e-14/1.083e4)
    
    bp = S.ObsBandpass('wfc3,ir,f098m')
    obs_hei = S.Observation(hei, bp)
    obs_cont = S.Observation(continuum, bp)
    
def compare_all():
    files=glob.glob('*asn.fits')
    files.sort()
    for file in files:
        print file
        mywfc3.bg.compare_filter_sensitivity(Rap=0.4, mAB=25, root=file)
        
    #    
    import pysynphot as S
    bp = {}
    width = {}
    for filt in ['F098M', 'F105W', 'F110W']:
        bp[filt] = S.ObsBandpass('wfc3,ir,%s' %(filt.lower()))
        ok = bp[filt].throughput > (0.5*bp[filt].throughput.max())
        width[filt] = np.max(bp[filt].wave[ok]) - np.min(bp[filt].wave[ok])
        
    for filt in ['F098M', 'F105W', 'F110W']:
        w_ratio = width[filt]/width['F098M']
        l_ratio = bp[filt].pivot() / bp['F098M'].pivot()
        print '%s %.3f %.3f %.3f' %(filt, w_ratio, l_ratio, np.sqrt(w_ratio)*l_ratio)
     
    os.system('head -1 ibl70k020_filter_SN.dat > SN_result.dat ; cat *SN.dat |grep -v "\#" >> SN_result.dat')
    x = catIO.Readfile('SN_result.dat')
    list = np.array(glob.glob('*SN.dat'))
    list.sort()
    
    ratio_f105w = x.obs_f105w / x.best_f098m
    bad = (x.obs_f105w / x.best_f105w) > 1 ## zodi underestimated
    ratio_f105w[bad] = (x.best_f105w / x.best_f098m)[bad]
    
    ratio_f110w = x.obs_f110w / x.best_f105w
    ratio_f110w[bad] = (x.best_f110w / x.best_f105w)[bad]
    
    optimal_f105w = np.median(x.best_f105w / x.best_f098m)
    optimal_f110w = np.median(x.best_f110w / x.best_f105w)
    
    ########### Make figure
    la, ra, ta, ba = 0.06, 0.06, 0.07, 0.13
    NX, NY = 1, 1
    aspect = (NY*(1+ta+ba))/((2.*NX)*(1+la+ra))
    dx = (1-la-ra)/2.    
    
    fig = unicorn.plotting.plot_init(xs=8, aspect=aspect, left=0.08, right=0.025, top=0.015, bottom=0.08, hspace=0.0, wspace=0.17, use_tex=True, NO_GUI=True)
    
    #### Demo
    ax = fig.add_subplot(121)
    mywfc3.bg.compare_filter_sensitivity(Rap=0.4, mAB=25, root='ibp329i[qs]', cax=ax)    
    
    #### Prefer F105W over F098M ? 
    ax2 = fig.add_subplot(122)
    ran, bins = (0.6,1.60), 50
    yh, xh, nn = ax2.hist(ratio_f105w[ratio_f105w >= 1], range=ran, bins=bins, alpha=0.5, color='green', normed=False, histtype='stepfilled')
    yh, xh, nn = ax2.hist(ratio_f105w[ratio_f105w < 1], range=ran, bins=bins, alpha=0.5, color='blue', normed=False, histtype='stepfilled')
    
    ax2.plot(optimal_f105w*np.ones(2), [0,42], color='green', linewidth=2, alpha=0.8)
    ax2.set_xlabel(r'S/N: F105W$_\mathrm{observed}$ / F098M$_\mathrm{No\,\,He}$')
    ax2.set_ylabel(r'$N_\mathrm{visit}$')
    
    yy, dy = 46, 4
    ax2.fill_between([ran[0],1.0], (yy+dy)*np.ones(2), (yy-dy)*np.ones(2), color='blue', edgecolor='None', linewidth=4, alpha=0.8)
    ax2.text(ran[0]+(1-ran[0])/2., yy, 'F098M', color='white', va='center', ha='center')
    ax2.fill_between([1.0,ran[1]], (yy+dy)*np.ones(2), (yy-dy)*np.ones(2), color='green', edgecolor='None', linewidth=4, alpha=0.8)
    ax2.text(1+(ran[1]-1)/2., yy, 'Prefer F105W', color='white', va='center', ha='center')

    ax2.set_xlim(ran)
    ax2.set_ylim((0,50))
    
    unicorn.plotting.savefig(fig, 'choose_Y_filter.pdf')
    
    #### Prefer F110W over F105W ? 
    # ax3 = fig.add_subplot(133)
    # yh, xh, nn = ax3.hist(ratio_f110w, range=ran, bins=bins, alpha=0.5, color='red', normed=False, histtype='stepfilled')
    # ax3.plot(optimal_f110w*np.ones(2), [0,42], color='red', linewidth=2, alpha=0.8)
    # ax3.set_xlabel(r'S/N: F110W$_\mathrm{observed}$ / F105W$_\mathrm{No\,\,He}$')
    # ax3.set_ylabel(r'$N_\mathrm{visit}$')
    # 
    # ax3.fill_between([ran[0],1.0], (yy+dy)*np.ones(2), (yy-dy)*np.ones(2), color='green', edgecolor='None', linewidth=4, alpha=0.8)
    # ax3.text(ran[0]+(1-ran[0])/2., yy, 'F105W', color='white', va='center', ha='center')
    # ax3.fill_between([1.0,ran[1]], (yy+dy)*np.ones(2), (yy-dy)*np.ones(2), color='red', edgecolor='None', linewidth=4, alpha=0.8)
    # ax3.text(1+(ran[1]-1)/2., yy, 'Prefer F110W', color='white', va='center', ha='center')
    # 
    # ax3.set_xlim(ran)
    # ax3.set_ylim((0,50))
    
    unicorn.plotting.savefig(fig, 'choose_Y_filter.pdf')
    
    yh, xh, nn = plt.hist(ratio_f110w, range=ran, bins=bins, alpha=0.5, color='red', normed=True, histtype='stepfilled')
    plt.plot(optimal_f110w*np.ones(2), [0,37], color='red', linewidth=2, alpha=0.8)
    #plt.hist(x.best_f105w / x.best_f098m, range=(0,2), bins=100, alpha=0.5, color='green')
    
    efficiency = (x.obs_f105w / x.best_f105w)**2
    
    plt.scatter(x.texp, x.obs_f105w / x.best_f098m, alpha=0.5)
    
def compare_filter_sensitivity(Rap=0.4, mAB=25, root='ibp329i[qs]', cax=None):
    import pysynphot as S
    import threedhst
    from threedhst import catIO
    
    zodi = S.FileSpectrum('/Users/brammer/WFC3/Backgrounds/Synphot/zodiacal_model_001.fits')
    nz = zodi.renorm(25, S.units.VegaMag, S.ObsBandpass("V"))
    bp = S.ObsBandpass('wfc3,ir,%s' %('F105W'.lower()))
    obs_zodi = S.Observation(nz, bp)
    zodi_refmag = 25+2.5*np.log10(2.84/30.63)
    zodi_level = 1 # e/s/pix
    refmag = 25+2.5*np.log10(obs_zodi.countrate()/zodi_level)
    nz = zodi.renorm(refmag, S.units.VegaMag, S.ObsBandpass("V"))
    
    ### Compare S/N in F105W / F098M
    thermal = 1.556/30.63
    dark = 1.532/30.63
    RN = 14.6
    
    #Rap = 0.4
    Area = np.pi*(Rap/0.128254)**2
    enclosed = 0.85
    
    #mAB = 25
    
    continuum = S.FlatSpectrum(mAB, fluxunits='ABMag')
    #continuum = S.FlatSpectrum(mAB+2.5*np.log10(3.713), fluxunits='STMag')
    #bp = S.ObsBandpass('wfc3,ir,%s' %(filter.lower()))
    
    obs = {}
    zodi = {}
    for band in ['f098m', 'f105w', 'f110w']:
        obs[band] = S.Observation(continuum, S.ObsBandpass('wfc3,ir,%s' %(band))).countrate()
        zodi[band] = S.Observation(nz, S.ObsBandpass('wfc3,ir,%s' %(band))).countrate()
        
    # orbit = catIO.Readfile('ibp329iqj_F105W_orbit.dat')
    # orbit = catIO.Readfile('ibp329isj_F105W_orbit.dat')
    # orbit.nexp = orbit.read*0+1
    
    #root='ibp329i[qs]'
    if 'asn' in root:
        outfile = root.replace('_asn.fits', '_filter_SN.pdf')
        asn = threedhst.utils.ASNFile(root)
        exposures = []
        for exp in asn.exposures:
            x = glob.glob('%s*orbit.dat' %(exp[:-1]))
            if len(x) > 0:
                exposures.append(x[0])
        #
        if len(exposures) == 0:
            return False
    else:
        exposures = glob.glob('%s*orbit.dat' %(root))
        outfile = 'filter_SN.pdf'
    
    os.system('cat %s > f105w_example.dat' %(' '.join(exposures)))
    orbit = catIO.Readfile('f105w_example.dat')
    
    #### Force F105W zodi = 1.
    orbit.bg = orbit.bg - orbit.zodi + 1.
    orbit.zodi = orbit.zodi*0.+1.
    
    ni = 0
    orbit.nexp = orbit.read*0
    for i in range(orbit.N):
        if orbit.read[i] == 1:
            ni += 1
        #
        orbit.nexp[i] = ni
    
    #
    t0 = 0
    time = orbit.seconds*1
    for i in range(1, orbit.N):
        if orbit.seconds[i] < orbit.seconds[i-1]:
            t0 += orbit.seconds[i-1]
            print t0
        #
        time[i] += t0
        
    orbit.read = np.arange(1,orbit.N+1)
    #dt = 100
    #time = orbit.read*dt
    #time = orbit.seconds
    dt = np.append(time[0], np.diff(time))
    
    source_signal = {}
    bkg_signal = {}
    he_bkg_signal = {}
    variance = {}
    ratio = {}
    sn_optimal = {}
    sn_observed = {}
    
    if cax is None:
        fig = unicorn.plotting.plot_init(xs=4, aspect=1, left=0.09, top=0.01, bottom=0.08, right=0.01, hspace=0.03, wspace=0., use_tex=True, NO_GUI=True)
        ax = fig.add_subplot(111)
    else:
        ax = cax
        
    colors = {'f098m':'blue', 'f105w':'green', 'f110w':'red'}
    alphas = {'f098m':0.5, 'f105w':0.8, 'f110w':0.5}
    max = 0

    ax.plot([0,10],[500,500], color='black', alpha=0.8, linestyle='--', linewidth=3, label='with observed He background')

    for band in ['f110w', 'f105w', 'f098m']:
        source_signal[band] = time*obs[band]*enclosed
        bkg_signal[band] = time*(zodi[band]*Area*orbit.zodi[0] + (thermal + dark)*Area)
        variance[band] = source_signal[band] + bkg_signal[band] + RN**2*Area*orbit.nexp
        #
        label = band.upper() + ', observed'*(band == 'f105w') + ', inferred'*(band != 'f105w')
        ax.plot(time, source_signal[band]/np.sqrt(variance[band]), linewidth=3, alpha=alphas[band], label=label, color=colors[band])
        sn_optimal[band] = (source_signal[band]/np.sqrt(variance[band]))[-1]
        sn_observed[band] = sn_optimal[band]
        max = np.maximum(max, sn_optimal[band])
        ratio[band] = 1
        #
        he_bkg_signal = bkg_signal[band]*1
        if 'f1' in band:
            he_bkg_signal += np.cumsum(dt*(orbit.bg-orbit.zodi)*Area)
            #
            he_variance = source_signal[band] + he_bkg_signal + RN**2*Area*orbit.nexp
            ax.plot(time, source_signal[band]/np.sqrt(he_variance), linewidth=3, alpha=alphas[band], color=colors[band], linestyle='--')
            #
            sn_observed[band] = (source_signal[band]/np.sqrt(he_variance))[-1]
            ratio[band] = sn_optimal[band] / sn_observed[band]
        #
        print '%s  %.0f %d S/N=%.1f eff=%.2f' %(band, time[-1], orbit.nexp[-1], sn_optimal[band], ratio[band]**2)
    #
    
    ax.legend(loc='upper left', prop={'size':9}, handlelength=2.5, frameon=False)
    ax.set_xlim(0,1.05*time.max())
    ax.set_ylim(0,1.2*max)
    ax.set_xlabel('time (s)')
    ax.set_ylabel(r'S/N($t$)')
    ax.text(0.95*time.max(),0.2*max,'%s\n' %(root.replace('_asn.fits','')) + r'$m_\mathrm{AB}=%d$, R=%.1f$^{\prime\prime}$' %(mAB, Rap), ha='right', va='center') #, bbox=dict(edgecolor='black', facecolor='white', pad=10))
    
    if cax is None:
        unicorn.plotting.savefig(fig, outfile)
        #
        fp = open(outfile.replace('.pdf','.dat') ,'w')
        fp.write('# texp Nexp best_f098m best_f105w obs_f105w best_f110w obs_f110w\n# %s\n' %(' '.join(exposures)))
        fp.write(' %d  %d  %.3f  %.3f %.3f  %.3f %.3f\n' %(time[-1], len(exposures), sn_optimal['f098m'], sn_optimal['f105w'], sn_observed['f105w'], sn_optimal['f110w'], sn_observed['f110w']))
        fp.close()
    
    return exposures, sn_optimal, sn_observed
    
    # zod = S.FileSpectrum('/Users/brammer/WFC3/Backgrounds/Synphot/zodiacal_model_001.fits')
    # nz = zod.renorm(SB, S.units.VegaMag, S.ObsBandpass("V"))
    # bp = S.ObsBandpass('wfc3,ir,%s' %(FILTER.lower()))
    # obs = S.Observation(nz, bp)
    
def raw2flt():
    """
    *** Don't have to bother, get this from the last read of the IMA file! ***
    
    Idea, just take the last read as the full integration and subtract out a
    combination of low-line and high-line sky images
    
    (still need to create optimal sky images but current could work for test)
    
    have error array include full poisson term?
    
    - let multidrizzle flag the CRs
    
    - ignore saturated for now...
    """
    raw = pyfits.open('ibt311j6q_raw.fits')
    raw = pyfits.open('icat19ueq_raw.fits') # barro
    raw = pyfits.open('ica531hgq_raw.fits') # glass
    flt = pyfits.open(raw.filename().replace('raw','flt'))
    
    GFLAT = pyfits.open(os.getenv('iref')+ raw[0].header['PFLTFILE'].split('$')[1])[1].data

    DFLAT = {}
    DFLAT['G102'] = pyfits.open(os.path.join(os.getenv('iref'), 'uc72113oi_pfl.fits'))[1].data  # F105W
    DFLAT['G141'] = pyfits.open(os.path.join(os.getenv('iref'), 'uc721143i_pfl.fits'))[1].data  # F140W
    FLAT = DFLAT[raw[0].header['FILTER']]
    
    NSAMP = raw[0].header['NSAMP']
    
    im_dark = pyfits.open(os.getenv('iref')+ raw[0].header['DARKFILE'].split('$')[1])
    dark = im_dark['SCI', 16-NSAMP+1].data-im_dark['SCI',16].data
    
    GAIN = raw[0].header['CCDGAIN']
    
    ### Subtract zeroth order read
    LAST = 4
    sci1 = raw['SCI', LAST].data - raw['SCI', NSAMP].data
    SAMPTIME = raw['SCI', LAST].header['SAMPTIME']
    
    ### Non-linearity correction
    im_lin = pyfits.open(os.getenv('iref')+ raw[0].header['NLINFILE'].split('$')[1])
    Fc = sci1
    for i in range(1,5):
        Fc += im_lin['COEF', i].data*sci1**i
    
    #### Imaging flat correction
    counts = ((Fc - dark)*GAIN)
    
    #### UNITCORR (e/s), FLATCORR
    new_flt = (counts/SAMPTIME/GFLAT)[5:-5,5:-5]
    new_err = np.sqrt(counts + 20**2)/SAMPTIME
    
    # sky = pyfits.open(os.getenv('THREEDHST')+'/CONF/sky.G141.set005.fits')[0].data
    # sky = pyfits.open(os.getenv('THREEDHST')+'/CONF/sky.G141.set002.fits')[0].data
    # sky = pyfits.open(os.getenv('THREEDHST')+'/CONF/sky_g102_f105w.fits')[0].data
    # 
    # flt_sub = pyfits.open('../PREP_FLT/'+raw.filename().replace('raw','flt'))
    
    #### Check IMA
    ima2 = pyfits.open('ibhj03xvq_ima.fits')
    flt = pyfits.open(ima.filename().replace('ima','flt'))
    
    NSAMP = ima[0].header['NSAMP']
    dq = ima2['dq',1].data*0
    for i in range(NSAMP):
        dq = dq | ima2['dq',i+1].data
        ds9.view(ima['dq',i+1].data)
     
def shadow_phase(fits='ib5x51l5q_flt.fits.gz', info=None, verbose=True):
    """
    Compute which half of an orbit will be in shadow for a particular 
    RA/Dec at a given time
    
    If you don't have a FITS file you can set fits=None and 
    info = (ra, dec, mjd)
    
    coe: 
    import astropy.coordinates as co
    import astropy.time
    import mywfc3.bg

    macs = co.ICRS('06h47m50.27s +70d14m55.0s')
    obs_date = astropy.time.Time('2014-05-16 22:56:58', scale='utc')
    obs_date = astropy.time.Time('2014-12-26 22:56:58', scale='utc')

    mywfc3.bg.shadow_phase(fits=None, info=(macs.ra.deg, macs.dec.deg, obs_date.mjd))
    
    ### Rodney
    import astropy.units as u
    
    a2744 = co.ICRS('3.588333333333E+00	-3.039725000000E+01', unit=(u.deg, u.deg))
    obs_date = astropy.time.Time('2014-06-19 22:56:58', scale='utc')
    mywfc3.bg.shadow_phase(fits=None, info=(a2744.ra.deg, a2744.dec.deg, obs_date.mjd))
    mywfc3.bg.shadow_phase(fits='ica9t3a4q_raw.fits', info=None)
    
    
    """
    import ephem
    import cPickle as pickle
    import astropy.coordinates as co
    import astropy.units as u
    
    if hasattr(co, 'ICRS'):
        icrs = co.ICRS
    else:
        icrs = co.ICRSCoordinates
    
    import subprocess
    #hjd, hra, hdec = np.loadtxt('/Users/brammer/WFC3/Backgrounds/Synphot/sun_coords.dat', unpack=True)
    # fp = open('/Users/brammer/WFC3/Backgrounds/Synphot/sun_coords.pkl','wb')
    # pickle.dump(hjd, fp)
    # pickle.dump(hra, fp)
    # pickle.dump(hdec, fp)
    # fp.close()
    
    fp = open('/Users/brammer/WFC3/Backgrounds/Synphot/sun_coords.pkl','rb')
    hjd = pickle.load(fp)
    hra = pickle.load(fp)
    hdec = pickle.load(fp)
    fp.close()
       
    # fits = 'ib5x51l5q_flt.fits.gz' # beginning
    # ra, ra_antisun = 53.2726125, 149.20238228, 233.160416666
    # fits = 'ib5x0blyq_flt.fits.gz' # end
    # ra, ra_antisun, ra_antiobj = 53.1604166667, 313.360343191, 233.160416666
    # fits = 'ib5x31mqq_flt.fits.gz' # middle, just end
    # ra, ra_antisun, ra_antiobj = 53.257475, 45.2557456737, 233.257475
    
    if fits is not None:
        if 'gz' in fits:
            p = subprocess.Popen('gunzip -c %s | dfits - | fitsort RA_TARG DEC_TARG EXPSTART | tail -1 ' %(fits), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            ra, dec, jd = np.cast[float](stdout.split())
        else:
            # head = pyfits.getheader(fits, ext=0)
            # ra, dec, jd = head['RA_TARG'], head['DEC_TARG'], head['EXPSTART']
            p = subprocess.Popen('dfits %s | fitsort RA_TARG DEC_TARG EXPSTART | tail -1 ' %(fits), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            #print stdout
            ra, dec, jd = np.cast[float](stdout.split()[1:])
            
    else:
       ra, dec, jd = info
       
    ra_sun = np.interp(jd, hjd-24.e5-0.5, hra)
    dec_sun = np.interp(jd, hjd-24.e5-0.5, hdec)
        
    eq_targ = icrs(ra=ra, dec=dec, unit=(u.deg, u.deg))
    eq_sun = icrs(ra=ra_sun, dec=dec_sun, unit=(u.deg, u.deg))
    
    # ra_antisun = ra_sun - 180
    # if ra_antisun < 0:
    #     ra_antisun = ra_sun + 180
    # 
    # ra_antiobj = ra - 180
    # if ra_antiobj < 0:
    #     ra_antiobj = ra + 180
            
    ra_antisun = (ra_sun + 180) % 360.
    ra_antiobj = (ra + 180) #% 360.
    
    #print ra, ra_antisun, ra_antiobj
    
    import mywfc3.zodi
    lat, lng = mywfc3.zodi.helio_lat_lng(ra, dec, jd+24.e5+0.5)
    
    #print ra, dec, ra_sun, dec_sun, lat, lng
    
    delta = (ra_antiobj - ra_antisun + 180) % 360 - 180
    
    #if (ra_antisun > ra) & (ra_antiobj > ra_antisun):
    if delta > 0:
        bright = 'Beginning'
    else:
        bright = 'End'
        
    
    if verbose:
        print 'File    Ra   Ra_antisun  Ra_antiobj   Bright  Delta  Ecl_Lng'
        print '%s   %.3f  %.3f  %.3f  %s  %6.1f %6.1f' %(fits, ra, ra_antisun, ra_antiobj, bright, delta, lng)
        
    return ra, ra_antisun, delta, lng, (bright == 'End')*1
        
    #return (bright == 'End')*1 #+ (np.abs(ra-ra_antisun) < 20)*2
    
    #return dec-sun_dec, ra-sun_ra
    
    #ra, dec = 34.440618, -5.1721396
    # eq = co.ICRSCoordinates(ra=ra, dec=dec, unit=(u.deg, u.deg))
    # equat = ephem.Equatorial(str(eq.ra.format(sep=':', unit=u.hour)), str(eq.dec.format(sep=':', unit=u.deg)), epoch=ephem.J2000)
    # eclip_obs = ephem.Ecliptic(equat)
    # 
    # eq = co.ICRSCoordinates(ra=ra_sun, dec=dec_sun, unit=(u.deg, u.deg))
    # equat = ephem.Equatorial(str(eq.ra.format(sep=':', unit=u.hour)), str(eq.dec.format(sep=':', unit=u.deg)), epoch=ephem.J2000)
    # eclip_sun = ephem.Ecliptic(equat)
    # 
    # return (eclip_obs.lat-eclip_sun.lat)/np.pi*180, (eclip_obs.lon-eclip_sun.lon)/np.pi*180
    
def test_phase():
    
    d = catIO.Readfile('master.dat')
    ha = np.loadtxt('master.sun_ha')
    names = []
    for n in d.name:
        names.append(n.split('j_')[0]+'q')
    
    names = np.array(names)
       
    uniq = np.unique(names)
    colors = ['red','blue', 'orange', 'green']
    for n in uniq:
        f = glob.glob(n+'_flt.fits*')[0]
        bright = mywfc3.bg.shadow_phase(fits=f, verbose=False)
        xx = names == n
        print n
        plt.scatter(ha[xx], (d.bg-d.zodi)[xx], color=colors[bright], alpha=0.1)
        
def check_ephem():
    """
    N. Grogin sent me the HST ephemeris for the first set of A2744 F105W 
    visits.  Check that my estimates of the shadow entry correspond to
    the entries in the ephemeris
    """
    import subprocess
    import astropy.time as t
    import astropy.units as u
    
    fits = 'ic8n04wpq_flt.fits.gz'
    
    files = glob.glob('ic8n07*flt.fits.gz')
    tstart = []
    tend = []
    for fits in files:
        if 'gz' in fits:
            p = subprocess.Popen('gunzip -c %s | dfits - | fitsort EXPSTART EXPEND | tail -1 ' %(fits), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            expstart, expend = np.cast[float](stdout.split())
        else:
            head = pyfits.getheader(fits, ext=0)
            expstart, expend = head['EXPSTART'], head['EXPEND']
        #
        tstart.append(expstart)
        tend.append(expend)
    
    tstart = t.Time(tstart, format='mjd', scale='utc')
    tend = t.Time(tend, format='mjd', scale='utc')
    #plt.plot_date(texp.plot_date, np.ones(2), color='blue', alpha=0.5, linestyle='-', linewidth=4)
    
    ### Read SHADOW ephemeris
    e_ttag, e_d, e_item, e_comment = np.loadtxt('shadow_ephem.dat', unpack=True, dtype=str)

    ttag = []
    for tt in e_ttag:
        ttag.append('2013:'+tt)
        
    tshadow = t.Time(ttag, scale='utc', format='yday')
    i = e_comment[0] == 'EXIT'
    tshad_entry = tshadow[i::2]
    tshad_exit = tshadow[i+1::2]
    
    dt = t.TimeDelta(120*u.second)
    
    ### Show exposures
    NEXP = len(tstart)
    for i in range(NEXP):
        texp = t.Time([tstart[i], tend[i]])
        plt.plot_date(texp.plot_date, np.ones(2)*10, color='blue', alpha=0.5, linestyle='-', linewidth=1)
        plt.fill_between(texp.plot_date, 0.5*np.ones(2), 1.5*np.ones(2), color='blue', alpha=0.5)
    
    plt.gcf().autofmt_xdate()
    
    ### Show SHADOW
    dt = t.TimeDelta(2*u.hour)
    ok = (tshad_entry > (tstart[0]-dt)) & (tshad_entry < (tend[-1]+dt))
    
    ix = np.arange(len(ok))[ok]
    for i in ix:
        tshad = t.Time([tshad_entry[i], tshad_exit[i]])
        #plt.plot_date(tshad.plot_date, np.ones(2)+1, color='black', alpha=0.5, linestyle='-', linewidth=1)
        plt.fill_between(tshad.plot_date, 1.6*np.ones(2), 2.1*np.ones(2), color='black', alpha=0.5)
        
    ### Set xlimits
    dt = t.TimeDelta(10*u.minute)
    xlim = t.Time([tstart[0]-dt, tend[-1]+dt])
    plt.xlim(xlim.plot_date)
    
    plt.ylim(0,2.6)
    
def future_ephem():
    """
    Compute shadow characterisics for upcoming visits
    """
    import subprocess
    import astropy.time as t
    import astropy.units as u
    import os
    import numpy as np
    
    os.system('cat ephem_upcoming_1.dat | sed "s/SAA /SAA/" | sed "s/EXT,L=/EXIT/" | sed "s/ENT,L=/ENTRY/" | sed "s/[\(\)]//g" |grep -v Slew | awk \'{print $1, $2, $3, $4}\' > ephem_upcoming_1.reform')

    os.system('cat ephem_upcoming_2.dat | grep -v "OCC" | sed "s/SAA /SAA/" | sed "s/EXT,L=/EXIT/" | sed "s/ENT,L=/ENTRY/" | sed "s/[\(\)]//g" |grep -v Slew | awk \'{print $1, $2, $3, $4}\' > ephem_upcoming_2.reform')
    
    os.system('cat ephem_past_1.dat | grep -v "OCC" | sed "s/SAA /SAA/" | sed "s/EXT,L=/EXIT/" | sed "s/ENT,L=/ENTRY/" | sed "s/[\(\)]//g" |grep -v Slew | awk \'{print $1, $2, $3, $4}\' > ephem_past_1.reform')

    os.system('cat ephem_upcoming_m0416.dat | grep -v "OCC" | sed "s/SAA /SAA/" | sed "s/EXT,L=/EXIT/" | sed "s/ENT,L=/ENTRY/" | sed "s/[\(\)]//g" |grep -v Slew | awk \'{print $1, $2, $3, $4}\' > ephem_upcoming_m0416.reform')

    os.system('cat ephem_jun17_13496.dat | grep -v "OCC" | sed "s/SAA /SAA/" | sed "s/EXT,L=/EXIT/" | sed "s/ENT,L=/ENTRY/" | sed "s/[\(\)]//g" |grep -v Slew | awk \'{print $1, $2, $3, $4}\' > ephem_jun17_13496.reform')

    os.system('cat ephem_jun17_13498.dat | grep -v "OCC" | sed "s/SAA /SAA/" | sed "s/EXT,L=/EXIT/" | sed "s/ENT,L=/ENTRY/" | sed "s/[\(\)]//g" |grep -v Slew | awk \'{print $1, $2, $3, $4}\' > ephem_jun17_13498.reform')

    os.system('cat ephem_upcoming_13504_aug14.dat | grep -v "OCC" | sed "s/SAA /SAA/" | sed "s/EXT,L=/EXIT/" | sed "s/ENT,L=/ENTRY/" | sed "s/[\(\)]//g" |grep -v Slew | awk \'{print $1, $2, $3, $4}\' > ephem_upcoming_13504_aug14.reform')
    
    e_ttag, e_d, e_item, e_comment = np.loadtxt('ephem_upcoming_13504_aug14.reform', unpack=True, dtype=np.str)
    
    #e_ttag, e_d, e_item, e_comment = np.loadtxt('ephem_past_1.reform', unpack=True, dtype=np.str)
    
    event = ['%s_%s' %(item, comment) for item, comment in zip(e_item, e_comment)]
    ttag = []
    for tt in e_ttag:
        ttag.append('2014:'+tt)
    
    # for tt in ttag:
    #     times = t.Time(tt, scale='utc', format='yday')
            
    times = t.Time(ttag, scale='utc', format='yday')
    
    skip = ['SAA27_EXIT', 'SAA27_ENTRY', 'SAA30_EXIT', 'TGT_AVD_ENTRY', 'FGS_AVD_EXIT', 'FGS_AVD_ENTRY', 'SHADOW_EXIT']
    keys = {'SAA30_ENTRY':0.2, 'TGT_AVD_EXIT':1, 'SHADOW_ENTRY':1.8}
    colors = {'SAA30_ENTRY':'green', 'TGT_AVD_EXIT':'blue', 'SHADOW_ENTRY':'black'}
    
    flt_start, flt_end = [], []
    
    #### Past exposures
    # flt_files = glob.glob('ic8n*flt*')
    # for file in flt_files:
    #     print file
    #     p = subprocess.Popen('gunzip -c %s | dfits - | fitsort EXPSTART EXPEND | tail -1 ' %(file), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #     stdout, stderr = p.communicate()
    #     expstart, expend = np.cast[float](stdout.split())
    #     flt_start.append(t.Time(expstart, scale='utc', format='mjd'))
    #     flt_end.append(t.Time(expend, scale='utc', format='mjd'))

    Nexp = len(flt_start)
    
    last = {}
    orbit_shadow = []
    t0 = times[0]
    for i in range(len(times)): #[:500]:
        if event[i] in skip:
            continue
        #
        ### Parse the duration
        dt_parse = np.cast[float](e_d[i][2:].split(':'))*np.array([u.day, u.hour, u.minute, u.second])
        dt = t.TimeDelta(dt_parse.sum())
        if event[i] == 'TGT_AVD_EXIT':
            acq = t.TimeDelta(8*u.minute)
        else:
            acq = t.TimeDelta(0*u.minute)
        #
        interval = t.Time([times[i]+acq, times[i]+dt])
        #
        plt.plot_date(interval.plot_date, np.ones(2)*keys[event[i]], color=colors[event[i]], alpha=0.5, linestyle='-', linewidth=1, marker='None')
        plt.fill_between(interval.plot_date, -0.5*np.ones(2)+keys[event[i]], 0.5*np.ones(2)+keys[event[i]], color=colors[event[i]], alpha=0.5)
        #
        last[event[i]] = interval
        if 'TGT' in event[i]:
            orbit_shadow.append((last['SHADOW_ENTRY'][1] - last['TGT_AVD_EXIT'][0]))
        #
        if (times[i]-t0) > 1*u.day:
            for j in range(Nexp):
                if (flt_start[j] > t0) & (flt_start[j] < times[i]):
                    exp_int = t.Time([flt_start[j], flt_end[j]])
                    plt.fill_between(exp_int.plot_date, -0.2*np.ones(2), 1.8*np.ones(2), color='red', alpha=0.5)
                    plt.gca().text(exp_int[1].plot_date, -0.3, flt_files[j].split('_flt')[0], rotation=90, size=8, ha='right', va='top')
            #
            plt.gcf().autofmt_xdate()
            plt.ylim(-1,3)
            t0.out_subfmt = 'date_hm'
            plt.title(t0.iso)
            t0.out_subfmt = 'date'
            plt.savefig('ephem_%s.png' %(t0.iso))
            print t0.iso
            plt.close()
            t0 = times[i]
    #
    plt.gcf().autofmt_xdate()
    plt.ylim(-1,3)
    t0.out_subfmt = 'date_hm'
    plt.title(t0.iso)
    t0.out_subfmt = 'date'
    plt.savefig('ephem_%s.png' %(t0.iso))
    plt.close()
    t0 = times[i]
    
    
    