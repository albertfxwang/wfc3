"""
Model grism spectra in individual FLTs

To run the fitting code you need to compile the cython code here with 

   ./setup_interp.py build_ext -i
   
"""
import os
import collections

import numpy as np
import scipy.ndimage as nd
import matplotlib.pyplot as plt

import astropy.io.fits as pyfits
from astropy.table import Table
import astropy.wcs as pywcs

import stwcs

### Helper functions from a document written by Pirzkal, Brammer & Ryan 
from .. import grism

photflam = {'F098M': 6.0501324882418389e-20, 'F105W': 3.038658152508547e-20, 'F110W': 1.5274130068787271e-20, 'F125W': 2.2483414275260141e-20, 'F140W': 1.4737154005353565e-20, 'F160W': 1.9275637653833683e-20}
photplam = {'F098M': 9864.722728110915, 'F105W': 10551.046906405772, 'F110W': 11534.45855553774, 'F125W': 12486.059785775655, 'F140W': 13922.907350356367, 'F160W': 15369.175708965562}
 
if False:
    import pysynphot as S
    n = 1.e-20
    spec = S.FlatSpectrum(n, fluxunits='flam')
    photflam = {}
    photplam = {}
    for filter in ['F098M', 'F105W', 'F110W', 'F125W', 'F140W', 'F160W']:
        bp = S.ObsBandpass('wfc3,ir,%s' %(filter.lower()))
        photplam[filter] = bp.pivot()
        obs = S.Observation(spec, bp)
        photflam[filter] = n/obs.countrate()
        
class GrismFLT(object):
    """
    Scripts for simple modeling of individual grism FLT images
    
    tbd: 
        o helper functions for extracting 2D spectra
        o lots of book-keeping for handling SExtractor objects & catalogs
        ...
        
    """
    def __init__(self, flt_file='ico205lwq_flt.fits', sci_ext=['SCI',1], refimage=None, segimage=None, refext=0, verbose=True, pad=100):
        """
        Todo: pull out stored PyFITS objects to improve pickling
        """
        ### Read the FLT FITS File
        self.flt_file = flt_file
        self.sci_ext = sci_ext
        
        #self.wcs = pywcs.WCS(self.im['SCI',1].header)
        #self.im = pyfits.open(self.flt_file)
        self.pad = pad
        self.read_flt()
        self.flt_wcs = stwcs.wcsutil.HSTWCS(self.im, ext=tuple(sci_ext))
        
        self.flt_wcs.naxis1 = self.im_header['NAXIS1']+2*self.pad
        self.flt_wcs.naxis2 = self.im_header['NAXIS2']+2*self.pad
        self.flt_wcs.wcs.crpix[0] += self.pad
        self.flt_wcs.sip.crpix[0] += self.pad
        self.flt_wcs.wcs.crpix[1] += self.pad
        self.flt_wcs.sip.crpix[1] += self.pad           
        
        #### Get PA
        self.get_pa()
        
        self.refimage = refimage
        self.refimage_im = None
        
        self.segimage = segimage
        self.segimage_im = None
        self.seg = np.zeros(self.im_data['SCI'].shape, dtype=np.float32)
        
        # print 'Break!'
        # return None
        
        if refimage is None:
            ### Use the FLT image itself
            self.filter = self.im[0].header['FILTER'].upper()
    
            ### Full science image scaled to F_lambda
            self.flam = self.im_data['SCI']*photflam[self.filter] #self.im_header['PHOTFLAM']
    
            ### xx Use full DQ mask for now
            self.dmask = self.im_data['DQ'] == 0
        else:
            ### Case where FLT is a grism exposure and reference direct image provided
            if verbose:
                print '%s / Blot reference image: %s' %(self.flt_file, refimage)
            
            self.refimage_im = pyfits.open(self.refimage)
            self.filter = self.refimage_im[refext].header['FILTER'].upper()
                            
            self.flam = self.get_blotted_reference(self.refimage_im, segmentation=False)*photflam[self.filter]
            self.dmask = np.ones(self.flam.shape, dtype=bool)
        
        if segimage is not None:
            if verbose:
                print '%s / Blot segmentation image: %s' %(self.flt_file, segimage)
            
            self.segimage_im = pyfits.open(self.segimage)
            self.process_segimage()
        
        self.pivot = photplam[self.filter]
                        
        # This needed for the C dispersing function
        self.clip = np.cast[np.double](self.flam*self.dmask)
        
        ### Read the configuration file.  
        ## xx generalize this to get the grism information from the FLT header
        self.grism = self.im_header0['FILTER'].upper()
        self.conf = grism.aXeConf(conf_file='/Users/brammer/3DHST/Spectra/Work/CONF/%s.test41.gbb.conf' %(self.grism))
        
        self.conf.get_beams()
        
        ### Beam arrays and sensitivities
        # self.dxlam = collections.OrderedDict()
        # self.nx = collections.OrderedDict()
        # self.sens = collections.OrderedDict()
        # for beam in self.conf.orders:
        #     if self.conf.orders[beam] > 0:
        #         self.dxlam[beam] = np.arange(self.conf.conf['BEAM%s' %(beam)][0], self.conf.conf['BEAM%s' %(beam)][1], dtype=int)
        #         self.nx[beam] = int(self.dxlam[beam].max()-self.dxlam[beam].min())+1
        #         self.sens[beam] = Table.read('%s/%s' %(os.path.dirname(self.conf.conf_file), self.conf.conf['SENSITIVITY_%s' %(beam)]))
        #         self.sens[beam].wave = np.cast[np.double](self.sens[beam]['WAVELENGTH'])
        #         self.sens[beam].sens = np.cast[np.double](self.sens[beam]['SENSITIVITY'])
        
        ### full_model is a flattened version of the FLT image
        self.modelf = np.zeros(self.sh_pad[0]*self.sh_pad[1])
        self.model = self.modelf.reshape(self.sh_pad)
        self.idx = np.arange(self.modelf.size).reshape(self.sh_pad)
        
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
    
    def read_flt(self):
        ### Read extensions of the file, adding padding if necessary
        self.im = pyfits.open(self.flt_file)
        self.im_header0 = self.im[0].header.copy()
        self.im_header = self.im[tuple(self.sci_ext)].header.copy()
        
        self.sh_flt = list(self.im[tuple(self.sci_ext)].data.shape)
        self.sh_pad = [x+2*self.pad for x in self.sh_flt]
        
        slx = slice(self.pad, self.pad+self.sh_flt[1])
        sly = slice(self.pad, self.pad+self.sh_flt[0])

        self.im_data = {}
        for ext in ['SCI', 'ERR', 'DQ']:
            iext = (ext, self.sci_ext[1])
            self.im_data[ext] = np.zeros(self.sh_pad, dtype=self.im[iext].data.dtype)
            self.im_data[ext][sly, slx] = self.im[iext].data*1
        
        self.im_data['DQ'] = self.unset_dq_bits(dq=self.im_data['DQ'], okbits=32+64+512)
        ### Negative pixels
        neg_pixels = self.im_data['SCI'] < -3*self.im_data['ERR']
        self.im_data['DQ'][neg_pixels] |= 1024
        self.im_data['SCI'][self.im_data['DQ'] > 0] = 0
        
        self.im_data_sci_background = False
        
    def clean_for_mp(self):
        """
        zero out io.fits objects to make suitable for multiprocessing parallelization
        """
        self.im = None
        self.refimage_im = None
        self.segimage_im = None
    
    def re_init(self):
        """
        Open io.fits objects again
        """
        self.im = pyfits.open(self.flt_file)
        if self.refimage:
            self.refimage_im = pyfits.open(self.refimage)
        if self.segimage:
            self.segimage_im = pyfits.open(self.segimage)
            
    def get_pa(self):
        h = self.im_header
        self.pa = 90 - np.arctan2(h['CD2_2'], h['CD2_1'])/np.pi*180
        
    def unset_dq_bits(self, dq=None, okbits=32+64+512, verbose=False):
        """
        Unset bit flags from a DQ array
        
        32, 64: these pixels seem OK
           512: blobs not relevant for grism exposures
        """
        bin_bits = np.binary_repr(okbits)
        n = len(bin_bits)
        for i in range(n):
            if bin_bits[-(i+1)] == '1':
                if verbose:
                    print 2**i
                
                dq -= (dq & 2**i)
        
        return dq
        
    def all_world2pix(self, ra, dec, idx=1, tolerance=1.e-4):
        """
        Handle awkward pywcs.all_world2pix for scalar arguments
        """
        if np.isscalar(ra):
            x, y = self.flt_wcs.all_world2pix([ra], [dec], idx, tolerance=tolerance, maxiter=100, quiet=True)
            return x[0], y[0]
        else:
            return self.flt_wcs.all_world2pix(ra, dec, idx, tolerance=tolerance, maxiter=100, quiet=True)
    
    def all_pix2world(self, x, y, idx=1):
        """
        Handle awkward pywcs.all_world2pix for scalar arguments
        """
        if np.isscalar(x):
            ra, dec = self.flt_wcs.all_pix2world([x], [y], idx)
            return ra[0], dec[0]
        else:
            return self.flt_wcs.all_pix2world(x, y, idx)
    
    def blot_catalog(self, catalog_table, ra='ra', dec='dec', sextractor=False):
        from astropy.table import Column
        
        if sextractor:
            ra, dec = 'X_WORLD', 'Y_WORLD'
            if ra.lower() in catalog_table.colnames:
                ra, dec = ra.lower(), dec.lower()
        
        
        tolerance=-4
        xy = None
        for wcs, wcsname in zip([self.flt_wcs, pywcs.WCS(self.im_header, relax=True)], ['HSTWCS', 'astropy.wcs']):
            if xy is not None:
                break
            for i in range(4):    
                try:
                    #xy = self.flt_wcs.all_world2pix(catalog_table[ra], catalog_table[dec], 1, tolerance=np.log10(tolerance+i))
                    xy = wcs.all_world2pix(catalog_table[ra], catalog_table[dec], 1, tolerance=np.log10(tolerance+i), quiet=True)
                    break
                except:
                    print '%s / all_world2pix failed to converge at tolerance = %d' %(wcsname, tolerance+i)
                
        sh = self.im_data['SCI'].shape
        keep = (xy[0] > 0) & (xy[0] < sh[1]) & (xy[1] > (self.pad-5)) & (xy[1] < (self.pad+self.sh_flt[0]+5))
        self.catalog = catalog_table[keep]
        
        for col in ['x_flt', 'y_flt']:
            if col in self.catalog.colnames:
                self.catalog.remove_column(col)
                
        self.catalog.add_column(Column(name='x_flt', data=xy[0][keep]))
        self.catalog.add_column(Column(name='y_flt', data=xy[1][keep]))
        
        if sextractor:
            self.catalog.rename_column('X_WORLD', 'ra')
            self.catalog.rename_column('Y_WORLD', 'dec')
            self.catalog.rename_column('NUMBER', 'id')
            
        return self.catalog
        
        if False:
            ### Compute full model
            self.modelf*=0
            for mag_lim, beams in zip([[10,24], [24,28]], ['ABCDEF', 'A']): 
                ok = (self.catalog['MAG_AUTO'] > mag_lim[0]) & (self.catalog['MAG_AUTO'] < mag_lim[1])
                so = np.argsort(self.catalog['MAG_AUTO'][ok])
                for i in range(ok.sum()):
                    ix = so[i]
                    print '%d id=%d mag=%.2f' %(i+1, self.catalog['NUMBER'][ok][ix], self.catalog['MAG_AUTO'][ok][ix])
                    for beam in beams:
                        self.compute_model(id=self.catalog['NUMBER'][ok][ix], x=self.catalog['x_flt'][ok][ix], y=self.catalog['y_flt'][ok][ix], beam=beam,  sh=[60, 60], verbose=True, in_place=True)
        
        #
        outm = self.compute_model(id=self.catalog['NUMBER'][ok][ix], x=self.catalog['x_flt'][ok][ix], y=self.catalog['y_flt'][ok][ix], beam='A',  sh=[60, 60], verbose=True, in_place=False)
        
    def show_catalog_positions(self, ds9=None):
        if not ds9:
            return False
        
        n = len(self.catalog)
        for i in range(n):
            ds9.set_region('circle %f %f 2' %(self.catalog['x_flt'][i], self.catalog['y_flt'][i]))
        
        if False:
            ok = catalog_table['MAG_AUTO'] < 23
            
    def process_segimage(self):
        """
        Blot the segmentation image
        """
        self.seg = self.get_blotted_reference(self.segimage_im, segmentation=True)
        
        self.seg_ids = np.cast[int](np.unique(self.seg)[1:])
        self.seg_flux = collections.OrderedDict()
        for id in self.seg_ids:
            self.seg_flux[id] = 0
            
        self.seg_mask = self.seg > 0
        for s, f in zip(self.seg[self.seg_mask], self.flam[self.seg_mask]):
            #print s
            self.seg_flux[int(s)] += f
        
        self.seg_flux = np.array(self.seg_flux.values())
    
    def get_segment_coords(self, id=1):
        """
        Get centroid of a given ID segment
        """
        yp, xp = np.indices(self.seg.shape)
        id_mask = self.seg == id
        norm = np.sum(self.flam[id_mask])
        xi = np.sum((xp*self.flam)[id_mask])/norm
        yi = np.sum((yp*self.flam)[id_mask])/norm
        return xi+1, yi+1, id_mask
    
    def align_bright_objects(self, flux_limit=4.e-18, xspec=None, yspec=None, ds9=None, max_shift=10, cutout_dimensions=[14,14]):
        ids = self.seg_ids[self.seg_flux > flux_limit]
        seg_dx = self.seg_ids*0.
        seg_dy = self.seg_ids*0.
        
        ccx0 = None
        for id in ids:
            xi, yi, id_mask = self.get_segment_coords(id)
            if (xi > 1014-201) | (yi < cutout_dimensions[0]) | (yi > 1014-20) | (xi < cutout_dimensions[1]):
                continue
                
            beam = BeamCutout(x=xi, y=yi, cutout_dimensions=cutout_dimensions, beam='A', conf=self.conf, GrismFLT=self) #direct_flam=self.flam, grism_flt=self.im)
            ix = self.seg_ids == id
            
            dx, dy, ccx, ccy = beam.align_spectrum(xspec=xspec, yspec=yspec)
            
            ### Try eazy template
            # sed = eazy.getEazySED(id-1, MAIN_OUTPUT_FILE='goodsn_3dhst.v4.1.dusty', OUTPUT_DIRECTORY='/Users/brammer/3DHST/Spectra/Release/v4.1/EazyRerun/OUTPUT/', CACHE_FILE='Same', scale_flambda=True, verbose=False, individual_templates=False)
            # dx, dy, ccx, ccy = beam.align_spectrum(xspec=sed[0], yspec=sed[1]/self.seg_flux[ix])
            plt.plot(ccx/ccx.max(), color='b', alpha=0.3)
            plt.plot(ccy/ccy.max(), color='g', alpha=0.3)
            
            if ccx0 is None:
                ccx0 = ccx*1.
                ccy0 = ccy*1.
                ncc = 0.
            else:
                ccx0 += ccx
                ccy0 += ccy
                ncc += 1
                
            if (np.abs(dx) > max_shift) | (np.abs(dy) > max_shift):
                print '[bad] ID:%d (%7.1f,%7.1f), offset=%7.2f %7.2f' %(id, xi, yi, dx, dy)
                #break
                continue
            
            seg_dx[ix] = dx
            seg_dy[ix] = dy
            print 'ID:%d (%7.1f,%7.1f), offset=%7.2f %7.2f' %(id, xi, yi, dx, dy)
            
            if ds9:
                beam.init_dispersion(xoff=0, yoff=0); beam.compute_model(beam.thumb, xspec=xspec, yspec=yspec); m0 = beam.model*1    
                beam.init_dispersion(xoff=-dx, yoff=-dy); beam.compute_model(beam.thumb, xspec=xspec, yspec=yspec); m1 = beam.model*1
                ds9.view(beam.cutout_sci-m1)
        
        ok = seg_dx != 0
        xsh, x_rms = np.mean(seg_dx[ok]), np.std(seg_dx[ok])
        ysh, y_rms = np.mean(seg_dy[ok]), np.std(seg_dy[ok])
        print 'dx = %7.3f (%7.3f), dy = %7.3f (%7.3f)' %(xsh, x_rms, ysh, y_rms)
        return xsh, ysh, x_rms, y_rms
    
    def update_wcs_with_shift(self, xsh=0.0, ysh=0.0, drizzle_ref=True):
        """
        Update WCS
        """
        import drizzlepac.updatehdr
        #drizzlepac.updatehdr.updatewcs_with_shift(self.im, self.im.filename(), rot=rot, scale=scale, xsh=xsh, ysh=ysh)
        #self.flt_wcs = stwcs.wcsutil.HSTWCS(self.im, ext=('SCI',1))
        
        if hasattr(self, 'xsh'):
            self.xsh += xsh
            self.ysh += ysh
        else:
            self.xsh = xsh
            self.ysh = ysh
            
        h = self.im_header
        rd = self.all_pix2world([h['CRPIX1'], h['CRPIX1']+xsh], [h['CRPIX2'], h['CRPIX2']+ysh])
        h['CRVAL1'] = rd[0][1]
        h['CRVAL2'] = rd[1][1]
        self.im[tuple(self.sci_ext)].header = h
        self.flt_wcs = stwcs.wcsutil.HSTWCS(self.im, ext=tuple(self.sci_ext))
        
        self.flt_wcs.naxis1 = self.im[sci_ext].header['NAXIS1']+2*self.pad
        self.flt_wcs.naxis2 = self.im[sci_ext].header['NAXIS2']+2*self.pad
        self.flt_wcs.wcs.crpix[0] += self.pad
        self.flt_wcs.sip.crpix[0] += self.pad
        self.flt_wcs.wcs.crpix[1] += self.pad
        self.flt_wcs.sip.crpix[1] += self.pad           
        
        if drizzle_ref:
            if (self.refimage) & (self.refimage_im is None):
                self.refimage_im = pyfits.open(self.refimage)
                 
            if self.refimage_im:
                self.flam = self.get_blotted_reference(self.refimage_im, segmentation=False)*photflam[self.filter]
                self.clip = np.cast[np.double](self.flam*self.dmask)
            
            if (self.segimage) & (self.segimage_im is None):
                self.segimage_im = pyfits.open(self.segimage)
            
            if self.segimage_im:
                self.process_segimage()
    
    def get_blotted_reference(self, refimage=None, segmentation=False, refext=0):
        """
        Use AstroDrizzle to blot reference / segmentation images to the FLT frame
        """
        import stwcs
        import astropy.wcs
        from drizzlepac import astrodrizzle
        
        #ref = pyfits.open(refimage)
        if refimage[refext].data.dtype != np.float32:
            refimage[refext].data = np.cast[np.float32](refimage[refext].data)
        
        refdata = refimage[refext].data
        if 'ORIENTAT' in refimage[refext].header.keys():
            refimage[refext].header.remove('ORIENTAT')
            
        if segmentation:
            ## todo: allow getting these from cached outputs for cases of very large mosaics            
            seg_ones = np.cast[np.float32](refdata > 0)-1
        
        # refdata = np.ones(refdata.shape, dtype=np.float32)
        # seg_ones = refdata
        
        #ref_wcs = astropy.wcs.WCS(refimage[refext].header)
        ref_wcs = stwcs.wcsutil.HSTWCS(refimage, ext=refext)
        #flt_wcs = stwcs.wcsutil.HSTWCS(self.im, ext=('SCI',1))
        flt_wcs = self.flt_wcs
        
        for wcs in [ref_wcs, flt_wcs]:
            if hasattr(wcs, 'idcscale'):
                if wcs.idcscale is None:
                    wcs.idcscale = np.sqrt(np.sum(wcs.wcs.cd[0,:]**2))*3600.
            else:
                wcs.idcscale = np.sqrt(np.sum(wcs.wcs.cd[0,:]**2))*3600.
            
            wcs.pscale = np.sqrt(wcs.wcs.cd[0,0]**2 + wcs.wcs.cd[1,0]**2)*3600.
            
            #print 'IDCSCALE: %.3f' %(wcs.idcscale)
            
        #print refimage.filename(), ref_wcs.idcscale, ref_wcs.wcs.cd, flt_wcs.idcscale, ref_wcs.orientat
            
        if segmentation:
            #print '\nSEGMENTATION\n\n',(seg_ones+1).dtype, refdata.dtype, ref_wcs, flt_wcs
            ### +1 here is a hack for some memory issues
            blotted_seg = astrodrizzle.ablot.do_blot(refdata+0, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=1, wcsmap=None)
            blotted_ones = astrodrizzle.ablot.do_blot(seg_ones+1, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=1, wcsmap=None)
            
            #blotted_seg = astrodrizzle.ablot.do_blot(refdata+0, ref_wcs, flt_wcs, 1, coeffs=True, interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None)
            #blotted_seg = astrodrizzle.ablot.do_blot(refdata, ref_wcs, flt_wcs, 1, coeffs=True, interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None)
            #blotted_ones = blotted_seg
            
            blotted_ones[blotted_ones == 0] = 1
            ratio = np.round(blotted_seg/blotted_ones)
            grow = nd.maximum_filter(ratio, size=3, mode='constant', cval=0)
            ratio[ratio == 0] = grow[ratio == 0]
            blotted = ratio
            
        else:
            #print '\nREFDATA\n\n', refdata.dtype, ref_wcs, flt_wcs
            blotted = astrodrizzle.ablot.do_blot(refdata, ref_wcs, flt_wcs, 1, coeffs=True, interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None)
        
        return blotted
        
    def compute_model(self, id=0, x=588.28, y=40.54, sh=[10,10], xspec=None, yspec=None, beam='A', verbose=False, in_place=True, outdata=None):
        """
        Compute a model spectrum, so simple!
        """
        
        from . import disperse
        import unicorn.utils_c
                 
        xc, yc = int(x), int(y)
        xcenter = x - xc
        
        ### Get dispersion parameters at the reference position
        dy, lam = self.conf.get_beam_trace(x=x-self.pad, y=y-self.pad, dx=self.conf.dxlam[beam]+xcenter, beam=beam)
        dyc = np.cast[int](dy+20)-20+1 # 20 for handling int of small negative numbers
        
        ### Account for pixel centering of the trace
        yfrac = dy-np.floor(dy)
        
        ### Interpolate the sensitivity curve on the wavelength grid. 
        ysens = lam*0
        so = np.argsort(lam)
        ysens[so] = unicorn.utils_c.interp_conserve_c(lam[so], self.conf.sens[beam]['WAVELENGTH'], self.conf.sens[beam]['SENSITIVITY'])/1.e17*np.abs(lam[1]-lam[0])
        if xspec is not None:
            yspec_int = ysens*0.
            yspec_int[so] = unicorn.utils_c.interp_conserve_c(lam[so], xspec, yspec)
            ysens *= yspec_int
        
        #print 'XXX', id, x, y, beam, ysens.max()
            
        x0 = np.array([yc, xc])
        slx = self.conf.dxlam[beam]+xc
        ok = (slx < self.sh_pad[1]) & (slx > 0)
        
        if in_place:
            #self.modelf *= 0
            outdata = self.modelf
        else:
            if outdata is None:
                outdata = self.modelf*0
        
        ### This is an array of indices for the spectral trace
        try:
            idxl = self.idx[dyc[ok]+yc,slx[ok]]
        except:
            if verbose:
                print 'Dispersed trace falls off the image: x=%.2f y=%.2f' %(x, y)
            
            return False
            
        ### Loop over pixels in the direct FLT and add them into a final 2D spectrum 
        ### (in the full (flattened) FLT frame)
        ## adds into the output array, initializing full array to zero would be very slow
        status = disperse.disperse_grism_object(self.clip, self.seg, id, idxl, yfrac[ok], ysens[ok], outdata, x0, np.array(self.clip.shape), np.array(sh), np.array(self.sh_pad))
                
        if not in_place:
            return outdata
        else:
            return True
    
    def fit_background(self, degree=3, sn_limit=0.1, pfit=None, apply=True, verbose=True, ds9=None):
        from astropy.modeling import models, fitting
        
        yp, xp = np.indices(self.sh_pad)
        xp  = (xp - self.sh_pad[1]/2.)/(self.sh_flt[1]/2)
        yp  = (yp - self.sh_pad[0]/2.)/(self.sh_flt[0]/2)
        
        if pfit is None:
            mask = (self.im_data['DQ'] == 0) & (self.model/self.im_data['ERR'] < sn_limit) & (self.im_data['SCI'] != 0) & (self.im_data['SCI'] > -4*self.im_data['ERR']) & (self.im_data['SCI'] < 6*self.im_data['ERR'])  
            poly = models.Polynomial2D(degree=degree)
            fit = fitting.LinearLSQFitter()
            pfit = fit(poly, xp[mask], yp[mask], self.im_data['SCI'][mask])
            pout = pfit(xp, yp)
            
            if ds9:
                ds9.view((self.im_data['SCI']-pout)*mask)
        else:
            pout = pfit(xp, yp)
            
        if apply:
            if self.pad > 0:
                slx = slice(self.pad, -self.pad)
                sly = slice(self.pad, -self.pad)

            else:
                slx = slice(0, self.sh_flt[1])
                sly = slice(0, self.sh_flt[0])
                
            self.im_data['SCI'][sly, slx] -= pout[sly, slx]
            self.im_data_sci_background = True
        
        if verbose:
            print 'fit_background, %s: p0_0=%7.4f' %(self.flt_file, pfit.parameters[0])
                
        self.fit_background_result = pfit #(pfit, xp, yp)
        
#### Testing
if False:
    import mywfc3.flt
    import threedhst.dq
    
    ds9 = threedhst.dq.myDS9()
    
    reload(mywfc3.flt.model); reload(mywfc3.flt); reload(mywfc3.grism)
    g = mywfc3.flt.model.GrismFLT(refimage='wfi2026-4536-05-235-f140w_drz_sci.fits')
    beam = mywfc3.flt.model.BeamCutout(conf=g.conf, cutout_dimensions=[20,20])
    flam_thumb = beam.get_flam_thumb(g.flam)
    
    beam.modelf[beam.idxl] += beam.ysens
    ds9.view(beam.modelf.reshape(beam.shg)) 
    
    beam.compute_model(flam_thumb)
    ds9.view(beam.modelf.reshape(beam.shg)) 
    
    gris = pyfits.open('ico205lxq_flt.fits')
    sly, slx = beam.get_slices()
    ds9.view(gris[1].data[sly, slx]-beam.modelf.reshape(beam.shg)) 
    
    ### 
    
class BeamCutout(object):
    """
    Cutout 2D spectrum from the full frame
    """
    def __init__(self, x=588.28, y=40.54, conf=None, cutout_dimensions=[10,10], beam='A', GrismFLT=None): #direct_flam=None, grism_flt=None, segm_flt=None):                 
        import unicorn.utils_c
        
        self.beam = beam
        self.x, self.y = x, y
        self.xc, self.yc = int(x), int(y)
        
        self.xcenter = self.x-self.xc
        
        if GrismFLT is not None:
            self.pad = GrismFLT.pad
        else:
            self.pad = 0
            
        self.dx = conf.dxlam[beam]
        
        self.cutout_dimensions = cutout_dimensions
        self.shd = np.array((2*cutout_dimensions[0], 2*cutout_dimensions[1]))
        self.lld = [self.yc-cutout_dimensions[0], self.xc-cutout_dimensions[1]]
        self.shg = np.array((2*cutout_dimensions[0], 2*cutout_dimensions[1] + conf.nx[beam]))
        self.llg = [self.yc-cutout_dimensions[0], self.xc-cutout_dimensions[1]+self.dx[0]]
        
        self.x_index = np.arange(self.shg[1])
        self.y_index = np.arange(self.shg[0])

        self.modelf = np.zeros(self.shg, dtype=np.double).flatten()
        self.model = self.modelf.reshape(self.shg)

        self.conf = conf
        self.beam = beam
        self.init_dispersion()
        
        self.thumb = None
        self.cutout_sci = None    
        self.shc = None
        self.cutout_seg = np.zeros(self.shg, dtype=np.float32)

        if GrismFLT is not None:
            self.thumb = self.get_flam_thumb(GrismFLT.flam)
            self.cutout_sci = self.get_cutout(GrismFLT.im_data['SCI'])*1
            self.cutout_dq = self.get_cutout(GrismFLT.im_data['DQ'])*1
            self.cutout_err = self.get_cutout(GrismFLT.im_data['ERR'])*1
            self.shc = self.cutout_sci.shape
            
            self.cutout_seg = self.get_cutout(GrismFLT.seg)*1
            
        # if direct_flam is not None:
        #     self.thumb = self.get_flam_thumb(direct_flam)
        # 
        # if grism_flt is not None:
        #     self.cutout_sci = self.get_cutout(grism_flt['SCI'].data)*1
        #     self.cutout_dq = self.get_cutout(grism_flt['DQ'].data)*1
        #     self.cutout_err = self.get_cutout(grism_flt['ERR'].data)*1
        #     self.shc = self.cutout_sci.shape
            #beam.cutout[beam.cutout_dq > 0] = 0
        
        #if segm_flt is not None:
         #   self.cutout_seg = self.get_cutout(segm_flt)*1
            
    def init_dispersion(self, xoff=0, yoff=0):
        """
        Allow for providing offsets to the dispersion function
        """
        import unicorn.utils_c
        
        self.dy, self.lam = self.conf.get_beam_trace(x=self.x-xoff-self.pad, y=self.y+yoff-self.pad, dx=self.conf.dxlam[self.beam]+self.xcenter-xoff, beam=self.beam)
        
        self.dy += yoff
        
        self.dyc = np.cast[int](self.dy)+1
        self.yfrac = self.dy-np.floor(self.dy)
        
        dl = np.abs(self.lam[1]-self.lam[0])
        self.ysens = unicorn.utils_c.interp_conserve_c(self.lam, self.conf.sens[self.beam]['WAVELENGTH'], self.conf.sens[self.beam]['SENSITIVITY'])/1.e17*dl
        
        self.idxl = np.arange(np.product(self.shg)).reshape(self.shg)[self.dyc+self.cutout_dimensions[0], self.dx-self.dx[0]+self.cutout_dimensions[1]]
        
    def get_flam_thumb(self, flam_full, xoff=0, yoff=0):
        dim = self.cutout_dimensions
        return np.cast[np.double](flam_full[self.yc+yoff-dim[0]:self.yc+yoff+dim[0], self.xc+xoff-dim[1]:self.xc+xoff+dim[1]])

    def compute_model(self, flam_thumb, id=0, yspec=None, xspec=None, in_place=True):
        from . import disperse
        import unicorn.utils_c
        
        x0 = np.array([self.cutout_dimensions[0], self.cutout_dimensions[0]])
        sh_thumb = np.array((self.shd[0]/2, self.shd[1]/2))  
        if in_place:
            self.modelf *= 0
            out = self.modelf
        else:
            out = self.modelf*0
            
        ynorm=1
        if xspec is not self.lam:
            if yspec is not None:
                ynorm = unicorn.utils_c.interp_conserve_c(self.lam, xspec, yspec)
        else:
            ynorm = yspec
            
        status = disperse.disperse_grism_object(flam_thumb, self.cutout_seg, id, self.idxl, self.yfrac, self.ysens*ynorm, out, x0, self.shd, sh_thumb, self.shg)
        
        if not in_place:
            return out
            
    def get_slices(self):
        sly = slice(self.llg[0], self.llg[0]+self.shg[0])
        slx = slice(self.llg[1], self.llg[1]+self.shg[1])
        return sly, slx
        
    def get_cutout(self, data):
        sly, slx = self.get_slices()
        return data[sly, slx]
        
    def make_wcs_header(self, data=None):
        import stwcs
        h = pyfits.Header()
        h['CRPIX1'] = self.cutout_dimensions[1]#+0.5
        h['CRPIX2'] = self.cutout_dimensions[0]#+0.5
        h['CRVAL1'] = self.lam[0]        
        h['CD1_1'] = self.lam[1]-self.lam[0]
        h['CD1_2'] = 0.
        
        h['CRVAL2'] = self.dy[0]
        h['CD2_2'] = 1.
        h['CD2_1'] = -(self.dy[1]-self.dy[0])
        
        if data is None:
            np.zeros(self.shg)
        
        data = hdul = pyfits.HDUList([pyfits.ImageHDU(data=data, header=h)])
        wcs = stwcs.wcsutil.HSTWCS(hdul, ext=0)
        wcs.pscale = np.sqrt(wcs.wcs.cd[0,0]**2 + wcs.wcs.cd[1,0]**2)*3600.
        
        return hdul[0], wcs
        
    def align_spectrum(self, xspec=None, yspec=None):
        """
        Try to compute alignment of the reference image using cross correlation
        """
        from astropy.modeling import models, fitting
        
        clean_cutout = self.cutout_sci*1.
        clean_cutout[self.cutout_dq > 0] = 0
        #max = np.percentile(clean_cutout[clean_cutout != 0], clip_percentile)
        #clean_cutout[(clean_cutout > max) | (clean_cutout < -3*self.cutout_err)] = 0
        clean_cutout[(clean_cutout < -3*self.cutout_err) | ~np.isfinite(self.cutout_err)] = 0.
        
        self.compute_model(self.thumb, xspec=xspec, yspec=yspec)
        
        ### Cross correlation
        cc = nd.correlate(self.model/self.model.sum(), clean_cutout/clean_cutout.sum())
        
        sh = cc.shape
        shx = sh[1]/2.; shy = sh[0]/2.

        yp, xp = np.indices(cc.shape)
        shx = sh[1]/2; shy = sh[0]/2
        xp = (xp-shx); yp = (yp-shy)

        cc[:,:shx-shy] = 0
        cc[:,shx+shy:] = 0
        ccy = cc.sum(axis=1)
        ccx = cc.sum(axis=0)
        
        #fit = fitting.LevMarLSQFitter()
        #mod = models.Polynomial1D(degree=6) #(1, 0, 1)
        fit = fitting.LinearLSQFitter()

        ix = np.argmax(ccx)
        p2 = models.Polynomial1D(degree=2)
        px = fit(p2, xp[0, ix-1:ix+2], ccx[ix-1:ix+2]/ccx.max())
        dx = -px.parameters[1]/(2*px.parameters[2])

        iy = np.argmax(ccy)
        py = fit(p2, yp[iy-1:iy+2, 0], ccy[iy-1:iy+2]/ccy.max())
        dy = -py.parameters[1]/(2*py.parameters[2])
        
        return dx, dy, ccx, ccy
        
                      