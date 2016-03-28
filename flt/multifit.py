"""
Fit spectra from multiple FLTs simultaneously
"""

import glob
import os
import time

import numpy as np
import matplotlib.pyplot as plt

import astropy.wcs
import astropy.io.fits as pyfits

from threedhst import catIO
import threedhst.dq

import mywfc3.flt
import mywfc3.flt.model
import mywfc3.flt.test_parallel

FLTs = None

def go():
    
    import os
    import mywfc3.flt.multifit
    reload(mywfc3.grism); reload(mywfc3.flt.model); reload(mywfc3.flt.test_parallel); reload(mywfc3.flt.multifit)
    
    ds9 = threedhst.dq.myDS9()
    
    os.chdir('/Users/brammer/3DHST/Spectra/Work/Erb/REDO_v2')
    files = ['ibr903dyq_flt.fits', 'ibr903e1q_flt.fits', 'ibr918e9q_flt.fits', 'ibr918ebq_flt.fits']
    self = mywfc3.flt.multifit.GroupFLT(files, refimage='Q0142-10-F140W_drz_sci.fits', segimage='Q0142-10-F140W_seg.fits', master_cat='Q0142-10-F140W.cat')
    
    os.chdir('/Users/brammer/3DHST/Spectra/Work/Oesch/TestFLT')
    files=glob.glob('*flt.fits')[:4]
    files.extend(glob.glob('*flt.fits')[12:16])
    
    files=glob.glob('*flt.fits')
    
    self = mywfc3.flt.multifit.GroupFLT(files, refimage='F160W_mosaic.fits', segimage='F160W_seg.fits', master_cat='/Users/brammer/3DHST/Spectra/Work/3DHST_Detection/GOODS-N_IR.cat')
    
    self.mp_compute_models(self.fit_data, initialize=True, verbose=1)
    print self.dt
    
    ### Testing
    id = 995
        
    dz = 0.003
    zr = [0.4, 3.5]
    
    cutout_dimensions=[10,10]

    beams = self.get_beams(id=id, cutout_dimensions=cutout_dimensions, ds9=ds9)
    mb = mywfc3.flt.multifit.MultiBeam(beams)
    
    refh, refwcs = mb.get_drizzle_wcs(center_wave=14000.0, dlam=46, NX=100, spatial_scale=1, NY=cutout_dimensions[0])
    refh['CRPIX2'] -= 2

    ref, refwcs = mb.beams[0].make_wcs_header(mb.beams[0].model)
    refh = ref.header
    #refh['CRVAL2'] *= 2# 1
    #refh['CD2_1'] = 0
    
    header, sci, wht, ctx = mb.drizzle(beam_attr='cutout_sci', pixfrac=0.3, ref_header=refh)
    header, contam, wht, ctx = mb.drizzle(beam_attr='contam', pixfrac=0.3, ref_header=refh)
    
    driz_beam = copy.copy(mb.beams[0])
    driz_beam.clean = (sci-contam).flatten()
    #driz_beam.cleanf = (sci-contam).flatten()
    driz_beam.ivar = wht
    mbd = mywfc3.flt.multifit.MultiBeam([driz_beam])
    
    zgrid, lnp, coeffs = mbd.fit_grid(zr=zr, dz=dz)
    plt.plot(zgrid, lnp-lnp.max())
        
    out = mb.fit_at_z(z=zgrid[np.argmax(lnp)], get_fit=True)
    mb.reshape_modelfit(out, attr='modelfit')
    #refh['CRPIX2'] += 0.1
    header, msci, mwht, mctx = mb.drizzle(beam_attr='modelfit', pixfrac=0.3, ref_header=refh)
    #plt.imshow(msci, interpolation='Nearest', vmin=-0.02, vmax=0.4)
        
def _loadFLT(file, refimage, segimage, cat):
    import mywfc3.flt.model
    
    if os.path.exists('%s.npy' %(file)):
        flt = np.load('%s.npy' %(file))[0]
    else:
        flt = mywfc3.flt.model.GrismFLT(flt_file=file, refimage=refimage, segimage=segimage)
        # xspec = np.arange(1.e4,1.8e4)
        # yspec = (xspec/1.6e4)**-0.4
        # xsh, ysh, x_rms, y_rms = flt.align_bright_objects(xspec=xspec, yspec=yspec, ds9=ds9)
        # flt.update_wcs_with_shift(xsh=xsh, ysh=ysh)
        c = flt.blot_catalog(cat, sextractor=False)
        
        ### Strip out io.fits objects
        flt.clean_for_mp()
    
    return flt
        
class GroupFLT():
    def __init__(self, files=[], refimage='Q0142-10-F140W_drz_sci.fits', segimage='Q0142-10-F140W_seg.fits', master_cat='Q0142-10-F140W.cat'):
        import multiprocessing as mp
        
        global FLTs
        
        self.FLTs = []
        self.files = files
        self.N = len(files)
        
        self.refimage = refimage
        self.segimage = segimage
        
        self.cat = catIO.Table(master_cat)
        self.cat.rename_column('X_WORLD', 'ra')
        self.cat.rename_column('Y_WORLD', 'dec')
        self.cat.rename_column('NUMBER', 'id')
        
        #### Read in the FLTs
        self._read_FLTs(self.files, initialize=True, parallel=(mp.cpu_count() > 1))
        
        # for i, file in enumerate(files):
        #     print '%d %s' %(i, file)
        #     if os.path.exists('%s.npy' %(file)):
        #         FLTs[i] = np.load('%s.npy' %(file))[0]
        #     else:
        #         flt = mywfc3.flt.model.GrismFLT(flt_file=file, refimage=refimage, segimage=segimage)
        #         # xspec = np.arange(1.e4,1.8e4)
        #         # yspec = (xspec/1.6e4)**-0.4
        #         # xsh, ysh, x_rms, y_rms = flt.align_bright_objects(xspec=xspec, yspec=yspec, ds9=ds9)
        #         # flt.update_wcs_with_shift(xsh=xsh, ysh=ysh)
        #         cat = flt.blot_catalog(self.cat, sextractor=False)
        #         self.FLTs.append(flt)
        # 
        # FLTs = self.FLTs
        
        self.init_fit_data()
    
    def _read_FLTs(self, files, initialize=True, parallel=True):
        import multiprocessing as mp
        global FLTs

        if initialize:
            FLTs = []
            
        if parallel:
            t0_pool = time.time()
            pool = mp.Pool(processes=mp.cpu_count())
    
            results = [pool.apply_async(_loadFLT, (self.files[i], self.refimage, self.segimage, self.cat)) for i in xrange(self.N)]
            pool.close()
            pool.join()
        
            for res in results:
                flt_i = res.get(timeout=1)
                FLTs.append(flt_i)
            
            t1_pool = time.time()
        else:
            for i in xrange(self.N):
                flt_i = _loadFLT(self.files[i], self.refimage, self.segimage, self.cat)
                FLTs.append(flt_i)
        
        ### Need to re-open the FLT files.  Don't bother with refimage/segimage 
        for flt in FLTs:
            flt.im = pyfits.open(flt.flt_file)
                
        self.FLTs = FLTs
            
    def init_fit_data(self, xspec=None, yspec=None, maglim={'A':25,'B':24,'C':22,'D':21,'E':21,'F':21}):
        if xspec is None:
            xspec = np.arange(9000,1.8e4,5)
        
        if yspec is None:
            yspec = (xspec/1.4e4)**-1 #-0.2
        
        keys = maglim.keys()
        keys.sort()
        
        self.fit_data = {}
        for i in range(len(self.cat)):
            beams = ''
            id_i = self.cat['id'][i] 
            mag_i = self.cat['MAG_AUTO'][i]
            for beam in keys:
                if mag_i <= maglim[beam]:
                    beams += beam
            
            if beams:
                self.fit_data[id_i] = {'xspec':xspec,'yspec':yspec,'beams':beams,'mag':mag_i}
            else:
                self.fit_data[id_i] = {'xspec':None,'yspec':None,'beams':beams,'mag':mag_i}
    
    def fix_fit_data(self, id):
        pass
    
    def get_beams(self, id=995, fcontam=1, cutout_dimensions=[14,14], ds9=None):
        
        beams = []
        for i, flt in enumerate(self.FLTs):
            if id not in flt.catalog['id']:
                continue
            
            ix = flt.catalog['id'] == id
            xi = flt.catalog['x_flt'][ix][0]
            yi = flt.catalog['y_flt'][ix][0]
        
            beam = mywfc3.flt.model.BeamCutout(x=xi, y=yi, cutout_dimensions=cutout_dimensions, beam='A', conf=flt.conf, direct_flam=flt.flam, grism_flt=flt.im)
            beam.cutout_seg = np.cast[np.float32](beam.get_flam_thumb(flt.seg))
            beam.id = id
            
            # import pysynphot as S
            # g = S.GaussianSource(1.e4, 1.4e4, 8)
            # beam.compute_model(beam.thumb, id=id, xspec=g.wave, yspec=g.flux)
            # beam.line_model = beam.model*1.
            
            beam.compute_model(beam.thumb, id=id, xspec=self.fit_data[id]['xspec'], yspec=self.fit_data[id]['yspec'])
            beam.contam = beam.get_cutout(flt.model) - beam.model
            
            beam.ivar = np.cast[np.float32](1/(beam.cutout_err**2+(fcontam*beam.contam)**2))
            beam.ivar[(beam.cutout_err == 0) | (beam.cutout_dq > 0) | (beam.cutout_err > 0.1)] = 0
            
            yp, xp = np.indices(beam.model.shape)
            #xp = (xp.flatten()-beam.cutout_dimensions[0])*(beam.lam[1]-beam.lam[0])+beam.lam[0]
            #xp = (xp.flatten()-beam.cutout_dimensions[0])*(beam.lam[1]-beam.lam[0])
            beam.xpf = (xp.flatten()-beam.cutout_dimensions[0]*1.)/100.
            beam.xp = beam.xpf.reshape(beam.shg)
            beam.cleanf = (beam.cutout_sci-beam.contam).flatten()
            beam.clean = beam.cleanf.reshape(beam.shg)
            
            beam.bg = 0.
            
            beams.append(beam)
            
            if ds9 is not None:
                ds9.frame(i+1)
                ds9.view((beam.cutout_sci-beam.contam)*(beam.cutout_dq == 0)) #, header=h.header)
        
        return beams
        
                
            # clf = sklearn.linear_model.LinearRegression()
            # status = clf.fit(A[ok,:], clean[ok]) #, sample_weight=ivar.flatten()[ok])
            # clf = sklearn.linear_model.Lasso()
            # status = clf.fit(A[ok,:], clean[ok], sample_weight=ivar.flatten()[ok])
        
    def fit_redshift(self, beams=[]):
        import sklearn.linear_model
        clf = sklearn.linear_model.Ridge(alpha=0.01)
        
        spec = np.loadtxt(os.getenv('THREEDHST')+'/templates/dobos11/SF0_0.emline.hiOIII.txt', unpack=True)
        
        A = None
        y = None
        ok = None
        ivar = None
        
        for beam in beams:
            beam.compute_model(beam.thumb, id=id)
            beam.cmodel = beam.modelf*1.
        
            oki = (beam.ivar.flatten() != 0) & (beam.cmodel > 0.05*beam.cmodel.max())
            
            if A is None:
                A = np.vstack([(beam.xpf*beam.cmodel), beam.cmodel, beam.cmodel*0., np.ones(beam.cmodel.size)]).T
                y = beam.clean*1
                ok = oki
                ivar = beam.ivar.flatten()
            else:
                Ap = np.vstack([(beam.xpf*beam.cmodel), beam.cmodel, beam.cmodel*0, np.ones(beam.cmodel.size)]).T
                A = np.append(A, Ap, axis=0)
                y = np.append(y, beam.clean*1)
                ok = np.append(ok, oki)
                ivar = np.append(ivar, beam.ivar.flatten())
                
        #### Start MCMC
        obj_fun = self._obj_zfit
        #obj_fun = _obj_zfit
        obj_args = [A, y, ivar, ok, beams, spec]

        NTHREADS=1
        ndim = 5

        NWALKERS, NSTEP = 100, 50
        #NWALKERS, NSTEP = 50, 50
        sampler = emcee.EnsembleSampler(NWALKERS, ndim, obj_fun, threads=NTHREADS, args=obj_args)
        
        z0 = np.random.rand(NWALKERS)*1.7+0.65
        z0 = 1.468+np.random.normal(size=NWALKERS)*0.1
        
        a0 = 0.25-np.random.rand(NWALKERS)*0.5
        b0 = 1.15-np.random.rand(NWALKERS)*0.3
        l0 = np.random.rand(NWALKERS)*2
        bg0 = np.random.normal(size=NWALKERS)*0.003
        
        p0 = np.array([z0, a0, b0, l0, bg0]).T
        result = sampler.run_mcmc(p0, NSTEP)
        
        param_names = ['z', 'a', 'b', 'line', 'bg']
        chain = unicorn.interlace_fit.emceeChain(chain=sampler.chain, param_names = param_names, sampler=sampler)
        
        xpix = (spec[0]*(1+chain.map[0]) - beam.lam[0])/(beam.lam[1]-beam.lam[0])/100.
        yspecf = (chain.map[1]*xpix + chain.map[2]) + chain.map[3]*spec[1]
        
        self.update_model(beam.id, xspec=spec[0]*(1+chain.map[0]), yspec=yspecf, verbose=True)
        for beam in beams:
            beam.compute_model(beam.thumb, id=beam.id, xspec=self.fit_data[id]['xspec'], yspec=self.fit_data[id]['yspec'])    
            beam.bg = chain.map[4]
        #
        self.show_beams(beams=beams, id=id, ds9=ds9, update_fit=False)
        
    @staticmethod
    def _obj_zfit(params, Ap, y, ivar, ok, beams, spec):
        
        z = params[0]
        if z < 0:
            return -np.inf
        
        coeffs = params[1:]
        
        for i, beam in enumerate(beams):
            beam.compute_model(beam.thumb, id=beam.id, xspec=spec[0]*(1+z), yspec=spec[1])
            if i == 0:
                lmodel = beam.modelf*1.
            else:
                lmodel = np.append(lmodel, beam.modelf*1.)
            
        A = Ap*1.
        A[:,2] = lmodel
        out = np.dot(A[ok,:], coeffs)
        lnp = -0.5*np.sum((out-y[ok])**2*ivar[ok])
        
        verb = '%.3f  ' %(z) + ' '.join(['%.1e' %c for c in coeffs]) + '%10.1f' %(lnp)
        print verb
        
        return lnp
        
    def update_model(self, id, xspec=None, yspec=None, verbose=False):
        """
        Update the FLT models with a new spectrum
        """
        fit_data_single = {}
        fit_data_single[id] = self.fit_data[id]
        fit_data_single[id]['yspec'] = fit_data_single[id]['yspec']*-1     
        self.mp_compute_models(fit_data=fit_data_single, initialize=False, verbose=verbose, parallel=False)
        
        fit_data_single[id]['xspec'] = xspec
        fit_data_single[id]['yspec'] = yspec
        self.mp_compute_models(fit_data=fit_data_single, initialize=False, verbose=verbose, parallel=False)
        
    def mp_compute_models(self, fit_data={}, initialize=True, verbose=True, parallel=True):
        import time
        import multiprocessing as mp
        
        if parallel:
            t0_pool = time.time()
            pool = mp.Pool(processes=mp.cpu_count())
        
            results=[pool.apply_async(_FLT_compute_model, (i, fit_data, initialize, verbose)) for i in xrange(self.N)]
            pool.close()
            pool.join()
        
            for res in results:
                i, modelf = res.get(timeout=1)
                self.FLTs[i].modelf = modelf
                self.FLTs[i].model = self.FLTs[i].modelf.reshape((1014,1014))
                #print FLTs[key].modelf.max()

            t1_pool = time.time()
        else:
            t0_pool = time.time()
            for i in xrange(self.N):
                out = _FLT_compute_model(i, fit_data, initialize, verbose)
            
            t1_pool = time.time()
                
        self.dt = t1_pool - t0_pool

class MultiBeam(object):
    def __init__(self, beams=[], model_min=0.05, ivar_min=[0.05,97]):
        self.beams = beams
        
        self.line_template_hi = np.loadtxt(os.getenv('THREEDHST')+'/templates/dobos11/SF0_0.emline.hiOIII.txt', unpack=True)
        
        #self.line_template_lo = np.loadtxt(os.getenv('THREEDHST')+'/templates/dobos11/SF0_0.emline.loOIII.txt', unpack=True)
        self.line_template_lo = np.loadtxt(os.getenv('THREEDHST')+'/templates/dobos11/SF0_0.emline.txt', unpack=True)
        
        #### Initialize multifit arrays
        A = None
        y = None
        mask = None
        ivar = None
        
        for beam in beams:
            beam.compute_model(beam.thumb, id=beam.id)
            beam.cmodel = beam.modelf*1.
        
            maski = (beam.ivar.flatten() != 0) & (beam.cmodel > model_min*beam.cmodel.max())
            
            if A is None:
                A = np.vstack([(beam.xpf*beam.cmodel), beam.cmodel, beam.cmodel*0., beam.cmodel*0., np.ones(beam.cmodel.size)]).T
                y = beam.clean*1
                mask = maski
                ivar = beam.ivar.flatten()
            else:
                Ap = np.vstack([(beam.xpf*beam.cmodel), beam.cmodel, beam.cmodel*0, beam.cmodel*0., np.ones(beam.cmodel.size)]).T
                A = np.append(A, Ap, axis=0)
                y = np.append(y, beam.clean*1)
                mask = np.append(mask, maski)
                ivar = np.append(ivar, beam.ivar.flatten())
        
        self.A = A
        self.y = y
        self.ivar = ivar
        self.mask = mask & (self.ivar > ivar_min[0]*np.percentile(self.ivar[mask], ivar_min[1]))
    
    def fit_at_z(self, z=0, get_fit=False):
        """
        TBD: add photometric constraints, at least from the direct image.
             - should help fits of bright continuum objects
             
        """
        import sklearn.linear_model
        
        for i, beam in enumerate(self.beams):
            lmodel_i = beam.compute_model(beam.thumb, id=beam.id, xspec=self.line_template_lo[0]*(1+z), yspec=self.line_template_lo[1], in_place=False)
            hmodel_i = beam.compute_model(beam.thumb, id=beam.id, xspec=self.line_template_hi[0]*(1+z), yspec=self.line_template_hi[1], in_place=False)
            if i == 0:
                lmodel = lmodel_i*1
                hmodel = hmodel_i*1
            else:
                lmodel = np.append(lmodel, lmodel_i)
                hmodel = np.append(hmodel, hmodel_i)
            
        Ap = self.A*1
        Ap[:,2] = lmodel
        Ap[:,3] = hmodel
        
        clf = sklearn.linear_model.LinearRegression(normalize=False, fit_intercept=False)
        status = clf.fit(Ap[self.mask,:], self.y[self.mask])
        
        #ls = np.linalg.lstsq(Ap[self.mask,:], self.y[self.mask])
        
        out = np.dot(Ap, clf.coef_)
        if get_fit:
            return out
            
        lnp = -0.5*np.sum((out-self.y)[self.mask]**2*self.ivar[self.mask])
        
        return lnp, clf.coef_
    
    def reshape_modelfit(self, out, attr='modelfit'):
        i0 = 0
        for beam in self.beams:
            i1 = beam.modelf.size
            beam.__setattr__(attr, out[i0:i0+i1].reshape(beam.shg))
            i0 += i1
            
        
    def fit_grid(self, zr=[0.7, 3.4], dz=0.01):
        """
        """
                
        zgrid = np.exp(np.arange(np.log(1+zr[0]), np.log(1+zr[1]), dz))-1
        
        
        nz = len(zgrid)
        lnp = np.zeros(nz)
        coeffs = np.zeros((nz, 5))
        for i in range(nz):
            lnp[i], coeffs[i,:] = self.fit_at_z(z=zgrid[i])
        
        return zgrid, lnp, coeffs
                    
    def get_drizzle_wcs(self, center_wave=1.4e4, dlam=40, NX=100, spatial_scale=1, NY=10):
        import astropy.io.fits as pyfits
        
        if False:
            center_wave=1.4e4; dlam=40; NX=100; spatial_scale=1; NY=10
        
        h = pyfits.ImageHDU(data=np.zeros((2*NY, 2*NX), dtype=np.float32))
        
        refh = h.header
        refh['CRPIX1'] = NX
        refh['CRPIX2'] = NY
        refh['CRVAL1'] = center_wave
        refh['CD1_1'] = dlam
        refh['CD1_2'] = 0.
        refh['CRVAL2'] = 0.
        refh['CD2_2'] = spatial_scale
        refh['CD2_1'] = 0.
        refh['RADESYS'] = ''
        
        ref_wcs = astropy.wcs.WCS(h.header)
        ref_wcs.pscale = np.sqrt(ref_wcs.wcs.cd[0,0]**2 + ref_wcs.wcs.cd[1,0]**2)*3600.
        
        return refh, ref_wcs
        
    def drizzle(self, beam_attr='clean', ref_header=None, pixfrac=1, ds9=None, update_fit=True, kernel='square'):
        from drizzlepac import astrodrizzle
              
        if ref_header is None:
            ref_header, ref_wcs = self.get_drizzle_wcs()
        else:
            ref_wcs = astropy.wcs.WCS(ref_header)
            ref_wcs.pscale = np.sqrt(ref_wcs.wcs.cd[0,0]**2 + ref_wcs.wcs.cd[1,0]**2)*3600.
                               
        sh = (ref_header['NAXIS2'], ref_header['NAXIS1'])
        
        outsci = np.zeros(sh, dtype=np.float32)
        outwht = np.zeros(sh, dtype=np.float32)
        outctx = np.zeros(sh, dtype=np.int32)
        
        for beam in self.beams:
            hdu, beam_wcs = beam.make_wcs_header(data=beam.clean)
            
            beam_data = beam.__getattribute__(beam_attr)
            if not (beam_data.dtype == np.float32):
                beam_data = np.cast[np.float32](beam_data) 
            
            if pixfrac == 0.:
                kernel = 'point'
                    
            astrodrizzle.adrizzle.do_driz(beam_data, beam_wcs, beam.ivar, ref_wcs, outsci, outwht, outctx, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=pixfrac, kernel=kernel, fillval=0, stepsize=10, wcsmap=None)   
            
            if ds9 is not None:
                ds9.view(outsci/ref_wcs.pscale**2)
        
        outsci /= ref_wcs.pscale**2
        outsci[outwht == 0] = 0
        return ref_header, outsci, outwht, outctx
        
        ##### optimal extractions
        wave = (np.arange(h['NAXIS1'])-h['CRPIX1'])*h['CD1_1'] + h['CRVAL1']
        spec = outsci.sum(axis=0)
        
        prof = outsci[:,(wave > 1.2e4) & (wave < 1.6e4)].sum(axis=1)
        prof /= prof.sum()
        
        num = prof*(outsci*outwht).T
        den = prof**2*outwht.T
        opt = num.T.sum(axis=0)/den.T.sum(axis=0)
        opt_var = 1./den.T.sum(axis=0)
        mnum = prof*(msci*outwht).T
        mopt = mnum.T.sum(axis=0)/den.T.sum(axis=0)
        
        #plt.plot(wave, spec, linestyle='steps-mid')
        plt.errorbar(wave, opt/driz_scale, np.sqrt(opt_var)/driz_scale, alpha=0.5, linestyle='steps-mid')
        #plt.plot(wave, opt/driz_scale, linestyle='steps-mid')
        plt.plot(wave, mopt/driz_scale, linestyle='steps-mid')
        plt.plot(wave, opt*0, linestyle='--', color='k')
        
        ### Test simple continuum + line fit where you can have a line that multiplies the continuum
        if update_fit:
            beam.compute_model(beam.thumb, id=id)
            beam.cmodel = beam.modelf*1.
        
            xspecf = np.arange(9000,1.8e4,1)
            yspecf = np.exp(-(xspecf-wave[np.argmax(opt)-1])**2/2./10**2)
            beam.compute_model(beam.thumb, id=id, xspec=xspecf, yspec=yspecf)
            beam.lmodel = beam.modelf*1.
            
            #yp, xp = np.indices(beam.model.shape)
            #xp = (xp.flatten()-beam.cutout_dimensions[0])*(beam.lam[1]-beam.lam[0])+beam.lam[0]
            #xp = (xp.flatten()-beam.cutout_dimensions[0])*(beam.lam[1]-beam.lam[0])
            ok = (beam.ivar.flatten() != 0) & (beam.cmodel > 0.0*beam.cmodel.max())
            A = np.vstack([(beam.xpf*beam.cmodel), beam.cmodel, beam.lmodel, np.ones(beam.lmodel.size)]).T
 
            #### Lstsq coeffs
            # out = np.linalg.lstsq(A[ok,:], clean[ok])
            # m = np.dot(A, out[0])

            #### Weighted Ridge coeffs            
            import sklearn.linear_model
            #clf = sklearn.linear_model.Ridge(alpha=0.01)
            #status = clf.fit(A[ok,:], beam.cleanf[ok], sample_weight=beam.ivar.flatten()[ok])
            clf = sklearn.linear_model.LinearRegression()
            status = clf.fit(A[ok,:], beam.cleanf[ok])
            m2 = np.dot(A, clf.coef_)
            
            xpix = (xspecf - beam.lam[0])/(beam.lam[1]-beam.lam[0])/100.
            yspecf = (clf.coef_[0]*xpix+clf.coef_[1]) + clf.coef_[2]*np.exp(-(xspecf-wave[np.argmax(opt)])**2/2./10**2)
            
            self.update_model(id, xspec=xspecf, yspec=yspecf, verbose=True)
                        
            for beam in beams:
                beam.compute_model(beam.thumb, id=id, xspec=self.fit_data[id]['xspec'], yspec=self.fit_data[id]['yspec'])    
                beam.bg = clf.coef_[3]
    
def _FLT_compute_model(ii, fit_data={}, initialize=True, verbose=True):
    import multiprocessing as mp
    #import multiprocessing.queues
    import time

    ### The only way I could figure out getting data in and out and passing through pool.apply_async
    global FLTs
    flt = FLTs[ii]

    ### Initialize
    if initialize:
        flt.modelf*=0
    
    keys = fit_data.keys()
    keys.sort()
    
    for i, id in enumerate(keys):
        if fit_data[id] is not None:
            if id in flt.catalog['id']:
                ix = flt.catalog['id'] == id
                for beam in fit_data[id]['beams']:
                    if verbose > 1:
                        print '%d %d %s id=%d mag=%.2f' %(ii, i+1, beam, id, fit_data[id]['mag'])
                    
                    flt.compute_model(id=id, x=flt.catalog['x_flt'][ix][0], y=flt.catalog['y_flt'][ix][0], beam=beam,  xspec=fit_data[id]['xspec'], yspec=fit_data[id]['yspec'], sh=[80, 80], verbose=False, in_place=True)
    
    if verbose:
        print '_FLT_compute_model - Done - %s' %(flt.flt_file)
                        
    # ### Loop through objects in given magnitude bins
    # #for mag_lim, beams in zip([[10,24], [24,28]], ['ABCDEF', 'A']): 
    # for mag_lim, beams in zip([[10,23]], ['A']): 
    #     ok = (flt.catalog['MAG_AUTO'] > mag_lim[0]) & (flt.catalog['MAG_AUTO'] < mag_lim[1])
    #     so = np.argsort(flt.catalog['MAG_AUTO'][ok])
    #     for i in range(ok.sum()):
    #         ix = so[i]
    #         if verbose:
    #             print '%d %s id=%d mag=%.2f' %(i+1, beam, flt.catalog['id'][ok][ix], flt.catalog['MAG_AUTO'][ok][ix])
    #         for beam in beams:
    #             flt.compute_model(id=flt.catalog['id'][ok][ix], x=flt.catalog['x_flt'][ok][ix], y=flt.catalog['y_flt'][ok][ix], beam=beam,  sh=[60, 60], verbose=False, in_place=True)

    #### For pool.apply_async or normal calls
    return (ii, flt.modelf)
    