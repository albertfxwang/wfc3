"""
Compute zodiacal light at a particular RA/Dec/Date

Get heliocentric lat/lon from lookup table
"""
import os

import numpy as np
import astropy.coordinates as co
import astropy.units as u
import pyfits

if hasattr(co, 'ICRS'):
    icrs = co.ICRS
else:
    icrs = co.ICRSCoordinates

def flt_zodi(image='ibhm04alq_flt.fits', verbose=True, pirzkal=False):
    from subprocess import Popen,PIPE

    if 'gz' in image:
        os.system('gunzip -c %s |dfits - | fitsort RA_TARG DEC_TARG EXPSTART DATE-OBS FILTER > /tmp/fits_head' %(image))
        params = open('/tmp/fits_head').readlines()[-1].split()
    else:
        stdout, stderr = Popen('gethead -x 0 %s RA_TARG DEC_TARG EXPSTART DATE-OBS FILTER' %(image), shell=True, stdout=PIPE).communicate()
        params = stdout.split()
    #
    ra, dec = np.cast[float](params[:2])
    jd = float(params[2]) + 2400000.5
    date = params[3]
    filter = params[4]
    
    if verbose:
        print ra, dec, date
    
    lat, lng = helio_lat_lng(ra, dec, jd)
    
    if pirzkal:
        return nor_zodi(lat, lng, filter=filter)
    else:
        return compute_zodi(ra, dec, jd, FILTER=filter, verbose=verbose)
    #
    # from subprocess import Popen,PIPE
    # stdout, stderr = Popen('dfits ibhm56f4q_raw.fits | fitsort RA_TARG DEC_TARG EXPSTART', shell=True, stdout=PIPE).communicate()
    # 
    # stdout, stderr = Popen('gethead -x 0 ibhm56f4q_raw.fits RA_TARG DEC_TARG EXPSTART', shell=True, stdout=PIPE).communicate()
    
    
def go():
    
    jd = 5.651374299768E+04 + 2400000.5
    ra, dec = 34.440618, -5.1721396

    print compute_zodi(ra, dec, jd, FILTER='F140W', verbose=True)
    
def compute_zodi(ra, dec, jd, FILTER='F140W', verbose=False):
    """
    Get the predicted Zodiacal surface brightness and then fold it through the Synphot WFC3 filters
    """
    import pysynphot as S
        
    thermal = {'F160W':0.134, 'F140W':0.069, 'F105W':0.051, 'F110W':0.05, 'F125W':0.052, 'G141':0.1, 'G102':0.04, 'F098M':0.05}
    
    lat, lng = helio_lat_lng(ra, dec, jd)
    SB = get_zodi_SB(lat, lng)
    if verbose:
        print 'Lat, Lng: %f, %f, SB=%.2f' %(lat, lng, SB+2.5*np.log10(0.128254**2))
        
    zod = S.FileSpectrum('/Users/brammer/WFC3/Backgrounds/Synphot/zodiacal_model_001.fits')
    nz = zod.renorm(SB, S.units.VegaMag, S.ObsBandpass("V"))
    bp = S.ObsBandpass('wfc3,ir,%s' %(FILTER.lower()))
    obs = S.Observation(nz, bp)
    
    return obs.countrate()+thermal[FILTER]
    
    
def helio_lat_lng(ra=0, dec=0, jd=0, fits=None):
    import ephem
    import cPickle as pickle
    #hjd, hra, hdec = np.loadtxt('/Users/brammer/WFC3/Backgrounds/Synphot/sun_coords.dat', unpack=True)
    fp = open('/Users/brammer/WFC3/Backgrounds/Synphot/sun_coords.pkl','rb')
    hjd = pickle.load(fp)
    hra = pickle.load(fp)
    hdec = pickle.load(fp)
    fp.close()
    
    ra_sun = np.interp(jd, hjd, hra)
    dec_sun = np.interp(jd, hjd, hdec)
    
    if fits is not None:
        head = pyfits.getheader(fits, ext=0)
        ra, dec, jd = head['RA_TARG'], head['DEC_TARG'], head['EXPSTART']
        
    #return dec-sun_dec, ra-sun_ra
    
    #ra, dec = 34.440618, -5.1721396
    eq = icrs(ra=ra, dec=dec, unit=(u.deg, u.deg))
    equat = ephem.Equatorial(str(eq.ra.format(sep=':', unit=u.hour)), str(eq.dec.format(sep=':', unit=u.deg)), epoch=ephem.J2000)
    eclip_obs = ephem.Ecliptic(equat)
    
    eq = icrs(ra=ra_sun, dec=dec_sun, unit=(u.deg, u.deg))
    equat = ephem.Equatorial(str(eq.ra.format(sep=':', unit=u.hour)), str(eq.dec.format(sep=':', unit=u.deg)), epoch=ephem.J2000)
    eclip_sun = ephem.Ecliptic(equat)
    
    #dlon = (eclip_obs.lon-eclip_sun.lon)/np.pi*180
    #print np.array([eclip_obs.lat, eclip_obs.lon, eclip_sun.lat, eclip_sun.lon])/np.pi*180
    
    dlon = ((eclip_obs.lon - eclip_sun.lon)/np.pi*180 + 180) % 360 - 180
    
    return (eclip_obs.lat-eclip_sun.lat)/np.pi*180, dlon
    
def get_zodi_SB(lat, lng):
    
    mat = np.loadtxt('/Users/brammer/WFC3/Backgrounds/Synphot/zodiac.txt')
    zlat = mat[0,1:]
    xlat = np.arange(len(zlat))
    zlng = mat[1:,0]
    xlng = np.arange(len(zlng))
    mat = mat[1:,1:]
    
    ilat = int(np.round(np.interp(np.abs(lat), zlat, xlat)))
    ilng = len(zlng)-1-int(np.round(np.interp(np.abs(lng), zlng[::-1], xlng)))
    
    val = mat[ilng, ilat]
    
    #### Vega V mag / pixel
    SB = 10-2.5*np.log10(val)-2.5*np.log10((1/3600.)**2)
    SB = 10-2.5*np.log10(val)-2.5*np.log10((0.128254/3600.)**2)
    
    return SB
#
def nor_zodi(la,lo,filter=None):
   """
   Nor's code to compute the Zodi background
   """
   import numpy as n
   la = la * 1.
   lo = lo * 1.
   mul = 1.
   if filter=="F098M":
       mul =  0.661509846941        
   elif filter=="F105W":
       mul = 1.22608999699
   elif filter=="F110W":
       mul = 2.05453276964  
   elif filter=="F125W":
       mul = 1.14614936917    
   elif filter=="F140W":
       mul = 1.46132279518  
   elif filter=="F160W":
       mul = 1.
   elif filter=="G141":
	mul = 2.08

   def c(la,lo):
        return n.cos(lo/180*n.pi) * n.cos(la/180*n.pi)
   def e(la,lo):
       return  n.arccos(c(lo, la))/n.pi*180
   def b(B):
       return 1.5 * (n.sqrt(1 + (B/1.5)**2) - 1 )

   a0 = 0.333885
   a1 = 0.0266097
   a2 = 0.628677
   a3 = 1.12545
   a4 = 1.70214
   a5 = 0.842976
   a6 = 0.00038706
   a7 = 2111.   

   res = a0 + a1 * (1 - n.cos(b(la)/180*n.pi)) + (a2 + a3 * c(lo, la) + a4 * (c(lo, la)**2) + a5* (c(lo, la)**3) ) * 10**(-n.sin(b(la)/180*n.pi)/(a6 * (e(lo, la) + a7))) 
   return res*mul
       
def compare_zodi():
    """
    Compare minimum background flux in a visit to the computed zodi
    """
    import mywfc3.zodi
    from threedhst import catIO
    import glob
    from mywfc3.utils import gzfile
    
    filter='G141'
    asns = glob.glob('*%s_orbit.png' %(filter))
    colors = np.array(['blue', 'red', 'orange'])
    #colors = np.array(['blue', 'white', 'orange'])
    # fig = unicorn.plotting.plot_init(xs=6, aspect=0.7, left=0.12, right=0.12)
    # ax = fig.add_subplot(111)
    bg_min = np.ones(len(asns))
    zodi_predicted = bg_min*1.
    nor = bg_min*1.
    for i, asn in enumerate(asns):
        root = asn.split('_')[0]
        print root
        os.system('cat %s*%s_orbit.dat > /tmp/%s.dat' %(root[:6], filter, root))
        bg1 = catIO.Readfile('/tmp/%s.dat' %(root), save_fits=False, force_lowercase=False)
        bg_min[i] = bg1.bg.min()
        files=glob.glob('%s*raw.fits*' %(root[:6]))
        zodi_predicted[i] = mywfc3.zodi.flt_zodi(gzfile(files[0]), verbose=False)
        nor[i] = mywfc3.zodi.flt_zodi(gzfile(files[0]), verbose=False, pirzkal=True)

    plt.scatter(bg_min, zodi_predicted, alpha=0.4, color='black', label='Synphot')
    plt.scatter(bg_min, nor, alpha=0.4, color='orange', label='Nor')

    plt.plot([0,4], [0,4], color='red', alpha=0.2, linewidth=3)
    plt.ylim(0.5,3.5); plt.xlim(0.5,3.5)
    plt.xlabel('Background, visit minimum')
    plt.ylabel('Predicted zodi')
    plt.legend(loc='lower right', prop={'size':10})
    
    plt.savefig('PredictedZodi_%s.pdf' %(filter))
    
    
    