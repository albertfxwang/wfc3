"""
Make a plot of the UVIS cosmic ray rate as a function of spacecraft 
suborbital position.  Shows that there is an enhancement in the CR
rate at longitudes near the magnetic poles, not just in the SAA.
"""
import numpy as np
from astropy.coordinates import Angle
import ephem

def archived_tle():
    """
    Return a list of HST TLEs from space-track
    https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2015-01-01--2016-12-31/NORAD_CAT_ID/20580/orderby/TLE_LINE1%20ASC/format/tle
    
    HST_TLE.txt: 
    https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2009-08-01--2016-12-31/NORAD_CAT_ID/20580/orderby/TLE_LINE1%20ASC/format/tle
    
    """
    from skyfield.api import load
    from skyfield.constants import AU_KM
    import astropy.time
    
    ts = load.timescale()
    eph = load('de421.bsp')
    earth = eph['earth']
    #sat = earth.satellite('\n'.join(lines))
    
    #lines = open('2015-16_TLEs.txt').readlines()
    lines = open('HST_TLE.txt').readlines()
    lines = [line.strip() for line in lines]

    N = len(lines)//2
    times = []
    strings = []
    
    for i in range(N):
        print(i,N)
        tle_str  = 'HST\n'+'\n'.join(lines[i*2:i*2+2])
        sat = earth.satellite(tle_str)
        t0 = astropy.time.Time(sat.epoch.utc_datetime())
        times.append(t0.mjd)
        strings.append(tle_str)
    
    return np.array(times), np.array(strings)
    
def go():
    """
    Read TLE file and save to a numpy file for faster I/O
    """
    tle_times, tle_strings = archived_tle()
    #np.save('2015-16_TLEs.npy', [tle_times, tle_strings])
    np.save('HST_TLE.npy', [tle_times, tle_strings])

def get_hst_positions(header, dt=10):
    """
    Get sublat/sublon of HST based on EXPSTART/STOP keywords in a 
    header
    """
    import ephem
    from astropy.coordinates import Angle
    import astropy.time
    import astropy.units as u
    
    #tle_times, tle_strings = np.load('2015-16_TLEs.npy')
    tle_times, tle_strings = np.load('HST_TLE.npy')

    tle_str = tle_strings[np.cast[float](tle_times) > header['EXPSTART']][0].split('\n')
    hst = ephem.readtle(tle_str[0], tle_str[1], tle_str[2])
    
    times = np.arange(header['EXPSTART'], header['EXPEND'], dt/86400.)
    pos = np.zeros((len(times), 2))
    for i, t in enumerate(times):
        tx = astropy.time.Time(t, format='mjd')
        hst.compute(tx.datetime)
        pos_i = np.array([Angle(hst.sublat*180/np.pi*u.deg).value, Angle(hst.sublong*180/np.pi*u.deg).wrap_at(360*u.deg).value])
        pos[i,:] = pos_i
        
    if pos[-1,1] < pos[0,1]:
        pos[pos[:,1] < pos[0,1], 1] += 360 
        
    return pos
        
    # t0x = astropy.time.Time(header['EXPSTART'], format='mjd')
    # hst.compute(t0x.datetime)
    # ll0 = np.array([Angle(hst.sublat*180/np.pi*u.deg).value, Angle(hst.sublong*180/np.pi*u.deg).wrap_at(360*u.deg).value])
    # 
    # t1x = astropy.time.Time(header['EXPEND'], format='mjd')
    # hst.compute(t1x.datetime)
    # ll1 = np.array([Angle(hst.sublat*180/np.pi*u.deg).value, Angle(hst.sublong*180/np.pi*u.deg).wrap_at(360*u.deg).value])
    # 
    # start_stop = np.array([ll0, ll1])
    # return start_stop
    
def test():
    jit = pyfits.open('/Users/brammer/Research/HST/UDS_Arc/LymanALpha/RAW/id8w04010_jit.fits')
    jit = pyfits.open('/Users/brammer/Research/HST/UDS_Arc/LymanALpha/RAW/id8w01010_jit.fits')
    jit = pyfits.open('/Users/brammer/Research/HST/UDS_Arc/LymanALpha/RAW/id8w02010_jit.fits')
    jit = pyfits.open('/Users/brammer/Research/HST/UDS_Arc/LymanALpha/RAW/id8w03010_jit.fits')

    jit = pyfits.open('/Users/brammer/3DHST/Spectra/Work/Grizli/WISP/RAW/id1kf0010_jit.fits')

    ix = 1

    flc = pyfits.open(glob.glob('{0}/{1}q_fl*.fits*'.format(os.path.dirname(jit.filename()), jit[ix].header['EXPNAME'][:-1]))[0])
    tab = astropy.table.Table.read(jit[ix])

    pl = plt.plot(tab['Longitude'], tab['Latitude'], linewidth=3, alpha=0.5)
    coo = get_hst_positions(flc[0].header)
    plt.plot(coo[:,1], coo[:,0], color=pl[0].get_color())

def show_darks():
    import costools
    import scipy.ndimage as nd
    import shapely
    from shapely.geometry import Polygon
    from descartes import PolygonPatch
    import astropy.io.fits as pyfits
    
    files = glob.glob("/grp/hst/wfc3s/bourque/blvs/ctecorr/*tmp.fits")
    files=glob.glob('*tmp.fits')
        
    ### Measure CR fraction
    files.sort()
    cr_fraction = np.zeros(len(files))
    headers = []
    
    for i, file in enumerate(files):
        im = pyfits.open(file)
        RN = 3.5
        dq = ((nd.minimum_filter(im['DQ',1].data, size=2) & 8192) > 0) & (im['SCI',1].data > 2*RN)
        cr_fraction[i] = dq.sum()*1./dq.size/im[0].header['EXPTIME']*1200
        headers.append(im[0].header)
        print(i, file, cr_fraction[i])
        
    plt.set_cmap('cubehelix')
    vm = [2, 5.4]
    import matplotlib.colors
    cm = plt.get_cmap()
    cnorm = matplotlib.colors.Normalize(vmin=vm[0], vmax=vm[1])
    
    fig = plt.figure(figsize=[8,4])
    import matplotlib.gridspec
    gs = matplotlib.gridspec.GridSpec(1,2, width_ratios=[1,2])
     
    ax = fig.add_subplot(gs[0])
    ax.hist(cr_fraction[:len(headers)]*100, range=[0, 7], bins=50, alpha=0.5, color='k')
    ax.set_xlabel('UVIS CR fraction\n% pixels in 1200 s')
    ax.set_ylabel('N')
    
    ax = fig.add_subplot(gs[1])
    for i, file in enumerate(files):
        print(i, file)
        #im = pyfits.open(file)
        #coo = get_hst_positions(im[0].header)
        coo = get_hst_positions(headers[i])
        
        N = coo.shape[0]//2
        
        # ax.plot(coo[:,1], coo[:,0], color='k', alpha=0.1)
        # ax.plot(coo[:,1]+360, coo[:,0], color='0.5', alpha=0.1)
        
        c_i = cm(cnorm(cr_fraction[i]*100))
        ax.plot(coo[:,1], coo[:,0], color=c_i, alpha=0.5, zorder=2)
        ax.plot(coo[:,1]+360, coo[:,0], color=c_i, alpha=0.5, zorder=2)
        
        #ax.scatter(coo[N:N+1,1], coo[N:N+1,0], c=[cr_fraction[i]*100], vmin=vm[0], vmax=vm[1], edgecolor='None', s=100, marker='s', alpha=0.8)
        #sc = ax.scatter(coo[N:N+1,1]+360, coo[N:N+1,0], c=[cr_fraction[i]*100], vmin=vm[0], vmax=vm[1], edgecolor='None', s=100, marker='s', alpha=0.8)
        
    for i in range(33):
        saa = np.array(costools.saamodel.saaModel(i))
        saa[:,1][saa[:,1] < 50] += 360
        #plt.plot(saa[:,1], saa[:,0], label=i, color='0.5', zorder=-1)
        
        poly = Polygon(saa[:,::-1])
        patch = PolygonPatch(poly, color='k', alpha=0.05, zorder=-1)
        ax.add_patch(patch)
        
    ax.set_xlim(180, 540)
    ax.set_ylim(-35, 35)
    ax.set_xlabel('SubLng')
    ax.set_ylabel('SubLat')
    ax.text(320, -25, 'SAA', ha='center', va='center', color='w')
    
    xarr = np.array(360+np.array([-150,-100,-50,0,50,100,150]))
    ax.set_xticks(xarr)
    ax.set_xticklabels(['150W', '100', '50', 'PM', '50E', '100', '150'])
    
    # magnetic poles
    ax.scatter([360-72.62, 360+107.38], [32, -32], marker='+', s=80, color='k')
    
    ax.grid()
    cb = plt.colorbar(sc)
    cb.set_label('CR fraction')
    
    gs.tight_layout(fig, pad=0.1)
        
