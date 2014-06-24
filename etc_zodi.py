"""
Copy code from PyETC to evaluate Zodi prediction as a specific time/coord
"""
from __future__ import division

import datetime
import os

import pyfits
import numpy as np

from pyetc.etc_engine.util import read_dict
from pyetc.etc_engine import exceptions
from pyetc.etc_web.etc import sky_angles as angles

#TODO: sky.dat needs to migrate to CDBS & have a versioned name
MODELS = read_dict(os.path.join(os.path.dirname(__file__),
                                'data/sky.dat'))
                                
#
class ZodiacalLight():
    """Zodiacal light is modeled by a spectrum file that can be
    normalized in a variety of ways, including based on the angle
    between the sun and the observed coordinates.
    """

    specfile = "spec(%s)" %(MODELS['zodi'])
    ZODI_LEVELS = {  #Normalized to Johnson V
        'None': None,
        'Low': 23.3,
        'Average': 22.7,
        'High': 22.1}

    #TODO: handle the posfile in the same way as the other data files
    #  except that some variable expansion has to be done, unlike for
    #  the other files for which that is handled by pysynphot
    #posfile = MODELS['zodipos']
    
    def __init__(self, flt='ibrh04qmq_flt.fits', level=None):
        import pyfits
        
        self.level = level
        if level is not None:
            self.ztype = 'ZodiStandard'
        else:
            self.im = pyfits.open(flt)
            self.ra = self.im[0].header['RA_TARG']
            self.dec = self.im[0].header['DEC_TARG']
            self.date = self.im[0].header['DATE-OBS']
            self.ztype = 'ZodiCoords'
        
        self._set_norm()
        self.make_synexpr()
        
    def eval_filter(self, filter=None, verbose=True):
        import pysynphot as S
        
        if filter is None:
            filter = self.im[0].header['FILTER']
        
        bp = S.ObsBandpass('wfc3,ir,%s' %(filter.lower()))
        #sp = S.FileSpectrum(self.expr)
        sp = S.FileSpectrum('$PYSYN_CDBS/etc/background/zodiacal_model_001.fits')
        rn = sp.renorm(self.norm - 2.5*np.log10(0.128254**2), 'Vegamag', S.ObsBandpass('johnson,v'))
        obs = S.Observation(rn, bp)
        self.obs = obs
        self.countrate = obs.countrate()
        
        ### Fudge factors to agree with ETC
        if filter == 'G102':
            self.countrate *= 1.3019
        
        if filter == 'G141':
            self.countrate *= 1.2261
            
        if verbose:
            print '%s, %s: %.3f e/s' %(self.descr, filter, self.countrate)
            
    def make_synexpr(self):
        if self.norm is None:
            self.expr = None
        else:
            if self.ztype in ('ZodiStandard', 'ZodiMag', 'ZodiCoords'):
                self.expr = "rn(%s,band(johnson,v),%f,vegamag)"%(self.specfile,self.norm)
            elif self.ztype in ('ZodiMult'):
                self.expr = "%f*%s"%(self.norm,self.specfile)

        
    def _set_norm(self):
        #There are many ways to specify how the Zodiacal Light spectrum should be normalized
        if self.ztype == 'ZodiStandard':
            level = self.level
            self.norm = self.ZODI_LEVELS[level]
            self.descr = level
            
        # elif self.ztype == 'ZodiMag':
        #     self.norm = self.fd[self.ztype]
        #     self.descr = 'User Defined (normalized to %f V)'%self.norm
        #     
        # elif self.ztype == 'ZodiMult':
        #     self.norm = self.fd[self.ztype]
        #     self.descr = 'User Defined (scaled by %f)'%self.norm

        elif self.ztype == 'ZodiCoords':
            #zra, zdec = self.ra, self.dec
            #ra, dec = angles.radec_to_decimal( zra, zdec )
            ra, dec = self.ra, self.dec
            
            #year, month, day = self.fd['ZodiDate'], self.fd['ZodiMonth'], self.fd['ZodiDay']
            year, month, day = np.cast[int](self.date.split('-'))
            
            date = datetime.date(year, month, day)
            helong, helat, sunangle = angles.radec_date_to_helonglat(ra, dec, date)
            self.descr = ('%6f %6f %s/%s/%s (%f deg.)' % (ra, dec, month, day, year, sunangle))
                
            try:
                self.norm = angles.zodi_lookup(helong, helat)
            except exceptions.RangeError, e:
                e.addinfo("These values are derived from the RA, Dec, and Date or Sunangle. You must either select different values for these quantities, or a different normalization method for the zodiacal light.")
                raise
                
        else:
            raise exceptions.BadEnumeration('Unknown normalization type %s' % self.ztype)
                                


