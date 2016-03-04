"""
Model grism spectra in individual FLTs
"""
import os
import collections

import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Table
import scipy.ndimage as nd

from .. import grism

class DirectFLT(object):
    def __init__(self, file='ico205lwq_flt.fits'):
        self.im = pyfits.open(file)
        self.filter = self.im[0].header['FILTER']

        self.flam = self.im['SCI'].data * self.im[0].header['PHOTFLAM']
        self.mask = self.im['DQ'].data == 0

        #msk2 = self.im['SCI'].data > 0.8*self.im['ERR'].data
        #self.clip = np.cast[np.double](self.flam*self.mask*msk2)
        self.clip = np.cast[np.double](self.flam * self.mask)
        
        self.conf = grism.aXeConf(conf_file='/Users/brammer/3DHST/Spectra/Work/CONF/G141.test41.gbb.conf')
        
        # self.sens = Table.read('%s/%s' %(os.path.dirname(self.conf.conf_file), self.conf.conf['SENSITIVITY_A']))
        # self.sens.wave = np.cast[np.double](self.sens['WAVELENGTH'])
        # self.sens.sens = np.cast[np.double](self.sens['SENSITIVITY'])
        
        self.dxlam = collections.OrderedDict()
        self.nx = collections.OrderedDict()
        self.sens = collections.OrderedDict()
        #Table.read('%s/%s' %(os.path.dirname(self.conf.conf_file), self.conf.conf['SENSITIVITY_A']))
        
        for beam in self.conf.orders:
            if self.conf.orders[beam] > 0:
                self.dxlam[beam] = np.arange(self.conf.conf['BEAM%s' %(beam)][0], self.conf.conf['BEAM%s' %(beam)][1], dtype=int)
                self.nx[beam] = int(self.dxlam[beam].max()-self.dxlam[beam].min())+1
                self.sens[beam] = Table.read('%s/%s' %(os.path.dirname(self.conf.conf_file), self.conf.conf['SENSITIVITY_%s' %(beam)]))
                self.sens[beam].wave = np.cast[np.double](self.sens[beam]['WAVELENGTH'])
                self.sens[beam].sens = np.cast[np.double](self.sens[beam]['SENSITIVITY'])
        
        self.full_model = np.zeros(1014**2)
        self.idx = np.arange(1014**2).reshape((1014,1014))
        
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
        Compute a model spectrum
        """
        
        from . import disperse
        import unicorn.utils_c
                 
        xc, yc = int(x), int(y)
        
        dy, lam = self.conf.get_beam_trace(x=x, y=y, dx=self.dxlam[beam], beam=beam)
        dyc = np.cast[int](dy)+1

        yfrac = dy-np.floor(dy)
        ysens = unicorn.utils_c.interp_conserve_c(lam, self.sens[beam].wave, self.sens[beam].sens)/1.e17*(lam[1]-lam[0])
        ysens *= (lam/1.e4)**beta

        x0 = np.array([xc, yc])
        slx = self.dxlam[beam]+xc
        ok = (slx < 1014) & (slx > 0)
        
        idxl = self.idx[dyc[ok]+yc,slx[ok]]
        
        status = disperse.disperse_grism_object(self.clip, idxl, yfrac[ok], ysens[ok], self.full_model, x0, np.array(sh))
        
        return True
        
        if False:
            ysensl = ysens*fl
            ysensh = ysens*fh
        
            twod = np.zeros((sh[0]*2+int(np.ceil(dy.max())), self.nx + 2*sh[1]))
        
            twod[sh[0]+dyc,sh[1]+self.dxlam-self.dxlam[0]] = fl*ysens
            twod[sh[0]+dyc-1,sh[1]+self.dxlam-self.dxlam[0]] = fh*ysens
        
            ix, iy = np.indices(twod.shape)
        
            twodc = nd.convolve(twod, thumb)
        
            # compare
            gris = pyfits.open('ico205lxq_flt.fits')
            sly = slice(yc-sh[0], yc-sh[0]+twod.shape[0])
            slx = slice(xc-sh[1]+self.dxlam[0], xc-sh[1]+self.dxlam[0]+2*sh[1]+self.nx)
            cutout = gris['SCI'].data[sly, slx]
            #cutout = 0
        
            #return twodc, cutout
        
            if True:
                twodf = twod.flatten()
                ### Try to speed it up
                nx, ny = twod.shape
                idx = np.arange(nx*ny).reshape(twod.shape)
                idxl = idx[sh[0]+dyc,sh[1]+self.dxlam-self.dxlam[0]]
                idxh = idx[sh[0]+dyc-1,sh[1]+self.dxlam-self.dxlam[0]]
            
                idx = np.arange(1014**2).reshape((1014,1014))
                idxl = idx[dyc+yc,self.dxlam+xc]
                idxh = idx[dyc+yc-1,self.dxlam+xc]
            
                idxb = np.append(idxl, idxh)
                ysensb = np.append(ysensl, ysensh)
            
                full = idx.flatten()*0.
                fluxf = self.flam.flatten()
                for k in range(10):
                    t0 = time.time()
                    for i in range(-sh[1], sh[1]):
                        for j in range(-sh[0], sh[0]):
                            full[idxl+j*1014+i] += ysensl*self.flam[yc+j, xc+i]*self.mask[yc+j, xc+i]/1.e-17
                            full[idxl+(j-1)*1014+i] += ysensh*self.flam[yc+j, xc+i]*self.mask[yc+j, xc+i]/1.e-17
                
                    print time.time()-t0
                
                #
                import mywfc3.flt.disperse
                flam = np.cast[np.double](self.flam*self.mask)
                msk2 = self.im['SCI'].data > 0.8*self.im['ERR'].data
                flam *= msk2
            
                x0 = np.array([xc, yc])
                x0 = np.array([507,507])
            
                full2 = idx.flatten()*0.
                yfrac = dy-np.floor(dy)
                status = mywfc3.flt.disperse.disperse_grism_object(flam, idxl, yfrac, ysens, full2, x0, np.array(sh))
            