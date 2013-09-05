"""
Get spacecraft pointing geometry from Jitter files.
"""
import os

import pyfits
import numpy as np
import matplotlib.pyplot as plt

import astropy.time
import astropy.units as u

import datetime

def go():
    import mywfc3.jitter
    
    mywfc3.jitter.show_orbit_limbangle(asn = ['ib3701050', 'ib3701060'])
    mywfc3.jitter.show_orbit_limbangle(asn = ['ib3702050', 'ib3702060'])
    for i in range(28):
        mywfc3.jitter.show_orbit_limbangle(asn = ['ib37%02d050' %(i+1), 'ib37%02d060' %(i+1)])
    
    asn_files = glob.glob('ibh*030_*asn.fits')
    for asn in asn_files:
        root = asn.split('_asn')[0]
        mywfc3.jitter.show_orbit_limbangle(asn = [root, root.replace('030', '040')])
        
    #mywfc3.jitter.show_orbit_limbangle(asn = ['ibhj20030', 'ibhj20040'])

def show_orbit_limbangle(asn = ['ib3701050', 'ib3701060']):
    
    import scipy.ndimage as nd
    
    os.chdir('/Users/brammer/WFC3/Jitter')
    
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
        flt = pyfits.open('../Backgrounds/G141/FLT/%s_flt.fits.gz' %(expname))
        # if FLAT_FILE != flt[0].header['PFLTFILE']:
        #     FLAT_IMAGE = pyfits.open(os.path.join(os.getenv('iref'), flt[0].header['PFLTFILE'].split('$')[1]))[1].data[5:-5, 5:-5]
        #     FLAT_FILE = flt[0].header['PFLTFILE']
        #
        flt[1].data /= FLAT_F140W
        xc, yc, NY = 707, 507, 100
        subim = flt[1].data[yc-NY:yc+NY, xc-NY:xc+NY]
        med_background = np.median(subim)
        ax3.plot(((pstr-tstr).sec + ext.data['Seconds'])/60., ext.data['LimbAng']*0+med_background, color=colors[i], linewidth=2)
        #
        if i > 3:
            ax_im = fig.add_axes((0.05+0.30/2*(i % 4)*1.05,0.09+0.28*2+0.01,0.30/2.,0.30))
            ax_im.imshow(nd.convolve(flt[1].data, np.ones((2,2))/4.)/med_background, vmin=0.6, vmax=1/0.6, interpolation='gaussian')
            ax_im.set_xticklabels([]); ax_im.set_yticklabels([])
            trace = np.median(flt[1].data/med_background, axis=0)
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
    
    plt.savefig('%s_orbit.png' %(targname))
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
                llcrnrlon=-180, urcrnrlon=-180+480, lat_ts=20, resolution='c', ax=ax)
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
    xpt,ypt = map(lon, lat)
    map.plot(xpt, ypt, alpha=0.5, **kwargs)
    map.scatter(xpt[0], ypt[0], alpha=0.5, **kwargs)
    #date = datetime.utcnow()
    #CS=m.nightshade(date, alpha=0.2, color='black')
    
        