"""
Model grism spectra in individual FLTs
"""
import os
import collections

import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Table
import scipy.ndimage as nd

### Helper functions from a document written by Pirzkal, Brammer & Ryan 
from .. import grism

class DirectFLT(object):
    """
    Scripts for simple modeling of individual grism FLT images
    
    tbd: 
        o helper functions for extracting 2D spectra
        o lots of book-keeping for handling SExtractor objects & catalogs
        ...
        
    """
    def __init__(self, file='ico205lwq_flt.fits'):
        
        ### Read the FLT FITS File
        self.im = pyfits.open(file)
        self.filter = self.im[0].header['FILTER']
        
        ### Full science image scaled to F_lambda
        self.flam = self.im['SCI'].data * self.im[0].header['PHOTFLAM']
        
        ### xx Use full DQ mask for now
        self.mask = self.im['DQ'].data == 0
        
        # This needed for the C dispersing function
        self.clip = np.cast[np.double](self.flam * self.mask)
        
        ### Read the configuration file.  
        ## xx generalize this to get the grism information from the FLT header
        self.conf = grism.aXeConf(conf_file='/Users/brammer/3DHST/Spectra/Work/CONF/G141.test41.gbb.conf')
        
        ### Beam arrays and sensitivities
        self.dxlam = collections.OrderedDict()
        self.nx = collections.OrderedDict()
        self.sens = collections.OrderedDict()
        for beam in self.conf.orders:
            if self.conf.orders[beam] > 0:
                self.dxlam[beam] = np.arange(self.conf.conf['BEAM%s' %(beam)][0], self.conf.conf['BEAM%s' %(beam)][1], dtype=int)
                self.nx[beam] = int(self.dxlam[beam].max()-self.dxlam[beam].min())+1
                self.sens[beam] = Table.read('%s/%s' %(os.path.dirname(self.conf.conf_file), self.conf.conf['SENSITIVITY_%s' %(beam)]))
                self.sens[beam].wave = np.cast[np.double](self.sens[beam]['WAVELENGTH'])
                self.sens[beam].sens = np.cast[np.double](self.sens[beam]['SENSITIVITY'])
        
        ### full_model is a flattened version of the FLT image
        self.full_model = np.zeros(1014**2)
        self.idx = np.arange(1014**2).reshape((1014,1014))
        
        ## Testing
        if True:
            pass
            #status = self.compute_model(x=281.13, y=856.8, sh=[20,20])
            # for xs in range(50,901,50):
            #     for ys in range(25,976,50):
            #         print xs, ys
            #         status = self.compute_model(x=xs, y=ys, sh=[25,25], beam='A')
            #         status = self.compute_model(x=xs, y=ys, sh=[25,25], beam='B')
            #         status = self.compute_model(x=xs, y=ys, sh=[25,25], beam='C')
            #         status = self.compute_model(x=xs, y=ys, sh=[25,25], beam='D')      
            #status = self.compute_model(x=507, y=507, sh=[100,100])
            
    def compute_model(self, x=588.28, y=40.54, sh=[10,10], beta=0, beam='A'):
        """
        Compute a model spectrum, so simple!
        
        TBD:
            o segmentation images
            o arbitrary input spectral models
            ...
        """
        
        from . import disperse
        import unicorn.utils_c
                 
        xc, yc = int(x), int(y)
        
        ### Get dispersion parameters at the reference position
        dy, lam = self.conf.get_beam_trace(x=x, y=y, dx=self.dxlam[beam], beam=beam)
        dyc = np.cast[int](dy)+1
        
        ### Account for pixel centering of the trace
        yfrac = dy-np.floor(dy)
        
        ### Interpolate the sensitivity curve on the wavelength grid. 
        ### xx still doesn't take an arbitrary spectrum as input but assumes F_lambda**beta
        ysens = unicorn.utils_c.interp_conserve_c(lam, self.sens[beam].wave, self.sens[beam].sens)/1.e17*(lam[1]-lam[0])
        ysens *= (lam/1.e4)**beta
        
        x0 = np.array([xc, yc])
        slx = self.dxlam[beam]+xc
        ok = (slx < 1014) & (slx > 0)
        
        ### This is an array of indices for the spectral trace
        idxl = self.idx[dyc[ok]+yc,slx[ok]]
        
        ### Loop over pixels in the direct FLT and add them into a final 2D spectrum (in the full (flattened) FLT frame)
        ## need better handling 
        status = disperse.disperse_grism_object(self.clip, idxl, yfrac[ok], ysens[ok], self.full_model, x0, np.array(sh))
        
        return True
        
            