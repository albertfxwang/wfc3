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

no_newline = '\x1b[1A\x1b[1M'

def go():
    
    import os
    import glob
    import mywfc3.flt.multifit
    reload(mywfc3.grism); reload(mywfc3.flt.model); reload(mywfc3.flt.test_parallel); reload(mywfc3.flt.multifit)
    import threedhst.dq
    
    ds9 = threedhst.dq.myDS9()
    
    #### Q0142
    # os.chdir('/Users/brammer/3DHST/Spectra/Work/Erb/REDO_v2')
    # files = ['ibr903dyq_flt.fits', 'ibr903e1q_flt.fits', 'ibr918e9q_flt.fits', 'ibr918ebq_flt.fits']
    # self = mywfc3.flt.multifit.GroupFLT(files, refimage='Q0142-10-F140W_drz_sci.fits', segimage='Q0142-10-F140W_seg.fits', master_cat='Q0142-10-F140W.cat')
    # 
    #### GOODS-N
    os.chdir('/Users/brammer/3DHST/Spectra/Work/Oesch/TestFLT')
    files=glob.glob('*flt.fits')[:4]
    files.extend(glob.glob('*flt.fits')[12:16])
    #files=glob.glob('*flt.fits')
    
    # self0 = mywfc3.flt.model.GrismFLT(flt_file=files[0], refimage='F160W_mosaic.fits', segimage='F160W_seg.fits', pad=0)
    # self100 = mywfc3.flt.model.GrismFLT(flt_file=files[0], refimage='F160W_mosaic.fits', segimage='F160W_seg.fits', pad=100)
    # refimage='F160W_mosaic.fits'; segimage='F160W_seg.fits'
    
    #self = mywfc3.flt.multifit.GroupFLT(files, refimage='F160W_mosaic.fits', segimage='F160W_seg.fits', master_cat='/Users/brammer/3DHST/Spectra/Work/3DHST_Detection/GOODS-N_IR.cat', pad=200)
    self = mywfc3.flt.multifit.GroupFLT(files, refimage='F160W_mosaic.fits', segimage='F160W_seg_blot.fits', master_cat='/Users/brammer/3DHST/Spectra/Work/3DHST_Detection/GOODS-N_IR.cat', pad=150)
    
    ### MACS0416
    info = catIO.Table('files.info')
    ok = ((info['FILTER'] == 'G141')) # & (info['PA_V3'] > 200)
    ok = ((info['FILTER'] == 'G102') | (info['FILTER'] == 'G141')) & (info['PA_V3'] > 200)
    files = [file.split('.gz')[0] for file in info['FILE'][ok]]
    
    files = list(info['FILE'][info['FILTER'] == 'G141'])
    files.extend(list(info['FILE'][info['FILTER'] == 'G102']))
    files = [file.split('.gz')[0] for file in files]
    
    self = mywfc3.flt.multifit.GroupFLT(files, refimage='MACS0416-F140W_drz_sci.fits', segimage='/Users/brammer/Research/HST/FrontierFields/DLV_Catalog/M0416/macs0416_detection_seg_grow.fits', master_cat='/Users/brammer/Research/HST/FrontierFields/DLV_Catalog/M0416/macs0416_detection.cat', pad=150)
    
    #### UDF
    # os.chdir('/Users/brammer/3DHST/Spectra/Work/UDF_FLT')
    # files=glob.glob('ibhj3*flt.fits')[:4]
    # files=glob.glob('icoi1ag*flt.fits')
    # 
    # f = mywfc3.flt.model.GrismFLT(flt_file='icoi1ag2q_flt.fits', refimage='/Users/brammer/3DHST/Ancillary/UDF/XDF/hlsp_xdf_hst_wfc3ir-60mas_hudf_f105w_v1_drz.fits', segimage='GOODS-S_IR.sub.seg.fits', pad=100)
    # 
    # #files.extend(glob.glob('*flt.fits')[12:16])
    # #files=glob.glob('*flt.fits')
    # self = mywfc3.flt.multifit.GroupFLT(files, refimage='/Users/brammer/3DHST/Ancillary/UDF/XDF/hlsp_xdf_hst_wfc3ir-60mas_hudf_f105w_v1_drz.fits', segimage='GOODS-S_IR.sub.seg.fits', master_cat='/Users/brammer/3DHST/Spectra/Work/3DHST_Detection/GOODS-S_IR.cat')
    # 
    # ### G102
    # self.init_fit_data(maglim={'A':25,'B':24,'C':22,'D':21,'E':21,'F':-100})
    
    self.mp_compute_models(self.fit_data, initialize=True, verbose=1, parallel=True)
    print self.dt
    
    ### Testing
    id = 995
        
    dz = 0.003
    zr = [0.4, 3.5]
    
    ###
    import unicorn
    c, z, f = unicorn.analysis.read_catalogs('goodsn')
    c, z, f = unicorn.analysis.read_catalogs('goodss')
    cat_zsp_ids = c['id'][(z['z_spec'] > 0.7) & (z['z_spec'] < 2.3)]
    ids = []
    z_spec = []
    for id in self.FLTs[0].catalog['id']:
        if id in cat_zsp_ids:
            ids.append(id)
            z_spec.append(z['z_spec'][z['id'] == id][0])
            
    cutout_dimensions=[10,10]
    
    full_beams = self.get_beams(id=id, cutout_dimensions=cutout_dimensions, ds9=None)
    driz_list = []

    for beams in [full_beams[:8]]:#, full_beams[8:]]:
        mb = mywfc3.flt.multifit.MultiBeam(beams)
    
        #refh, refwcs = mb.get_drizzle_wcs(center_wave=14000.0, dlam=46, NX=100, spatial_scale=1, NY=cutout_dimensions[0])
        #refh['CRPIX2'] -= 2

        ref, refwcs = mb.beams[0].make_wcs_header(mb.beams[0].model)
        refh = ref.header
        #refh['CRVAL2'] *= 2# 1
        #refh['CD2_1'] = 0
    
        header, sci, wht, ctx = mb.drizzle(beam_attr='cutout_sci', pixfrac=0.3, ref_header=refh, ds9=None)
        header, contam, wht, ctx = mb.drizzle(beam_attr='contam', pixfrac=0.3, ref_header=refh, ds9=None)
        ds9.view(sci-contam)
    
        driz_beam = copy.copy(mb.beams[0])
        driz_beam.clean = np.roll(sci-contam, -1, axis=1)
        driz_beam.cleanf = driz_beam.clean.flatten()
        driz_beam.thumb *= 0.
        for ib, beam in enumerate(beams):
            driz_beam.thumb += mb.beams[ib].thumb/len(beams)
        
        #driz_beam.cleanf = (sci-contam).flatten()
        driz_beam.ivar = np.roll(wht, -1, axis=1).flatten()
        driz_list.append(driz_beam)
        
    mbd = mywfc3.flt.multifit.MultiBeam(driz_list)
    
    zgrid, lnp, coeffs = mbd.fit_grid(zr=zr, dz=dz)
    z_gris, lnp_max = get_parabola_max(zgrid, lnp)
    plt.plot(zgrid, lnp-lnp_max)
    plt.scatter(z_gris, 0, color='k')
    
    print 'dz/(1+z): %.4f' %((z_gris-zsp)/(1+zsp))
    plt.plot(zsp*np.ones(2), [(lnp-lnp.max()).min(), 0], color='k')
    
    out, Ap, coeff = mb.fit_at_z(z=z_gris, get_fit=True)
    xfit, yfit = mb.get_1d_spec(coeff, z_gris)
    #self.update_model(id, xspec=xfit, yspec=yfit)
    
    mb.reshape_modelfit(out, attr='modelfit')
    #refh['CRPIX2'] += 0.1
    header, msci, mwht, mctx = mb.drizzle(beam_attr='modelfit', pixfrac=0.3, ref_header=refh)
    #plt.imshow(msci, interpolation='Nearest', vmin=-0.02, vmax=0.4)

def maxwcs_cutout(files, refimage=None, output='subimg.fits', pad=200):
    """
    Create a cutout image that encompases a list of input images
    """
    from drizzlepac import astrodrizzle
    import stwcs
    #self = mywfc3.flt.multifit.GroupFLT(files, refimage='F160W_mosaic.fits', segimage='/Users/brammer/3DHST/Spectra/Work/3DHST_Detection/goodsn_3dhst.v4.0.F160W_seg.fits', master_cat='/Users/brammer/3DHST/Spectra/Work/3DHST_Detection/GOODS-N_IR.cat', pad=200
    
    import astropy.wcs as pywcs
    list_of_wcsobj = []
    for file in files:
        im = pyfits.open(file)
        list_of_wcsobj.append(pywcs.WCS(im['SCI',1].header))
        
    imref = pyfits.open(refimage)
    ref_wcs = pywcs.WCS(imref[0].header)
        
    out_wcs = stwcs.distortion.utils.output_wcs(list_of_wcsobj, ref_wcs=None, owcs=None, undistort=True)
    
    for wcs in [ref_wcs, out_wcs]:
        wcs.pscale = np.sqrt(wcs.wcs.cd[0,0]**2 + wcs.wcs.cd[1,0]**2)*3600.

    map = astrodrizzle.ablot.wcs_functions.WCSMap(ref_wcs, out_wcs)
    px = np.array([0,1,1,0])*out_wcs._naxis1
    py = np.array([0,0,1,1])*out_wcs._naxis2

    px += (px > 0)*pad - (px == 0)*pad
    py += (py > 0)*pad - (py == 0)*pad
    ppx, ppy = np.cast[int](map.backward(px, py))
    
    slx = slice(ppx.min(), ppx.max())
    sly = slice(ppy.min(), ppy.max())
    
    subim = imref[0].data[sly, slx]
    sub_wcs = ref_wcs.slice((sly, slx))
    subh = sub_wcs.to_header()
    
    for i in range(2):
        for j in range(2):
            pc = 'PC%d_%d' %(i+1, j+1)
            cd = 'CD%d_%d' %(i+1, j+1)
            if pc in subh:
                subh[cd] = subh[pc]
                subh.remove(pc)
            else:
                subh[cd] = 0
    
    pyfits.writeto(output, data=subim, header=subh, clobber=True)
        
def get_parabola_max(x, y):
    from astropy.modeling import models, fitting
    
    fit = fitting.LinearLSQFitter()

    ix = np.argmax(y)
    if ix == 0:
        ix = 1
    
    if ix == (len(y)-1):
        ix -= 1
        
    p2 = models.Polynomial1D(degree=2)
    px = fit(p2, x[ix-1:ix+2], y[ix-1:ix+2]/y.max())
    dx = -px.parameters[1]/(2*px.parameters[2])
    return dx, px(dx)*y.max()
    
    
def _loadFLT(file, refimage, segimage, pad, cat):
    import mywfc3.flt.model
    
    if os.path.exists('%s.npy' %(file)):
        flt = np.load('%s.npy' %(file))[0]
    else:
        flt = mywfc3.flt.model.GrismFLT(flt_file=file, refimage=refimage, segimage=segimage, pad=pad)
        # xspec = np.arange(1.e4,1.8e4)
        # yspec = (xspec/1.6e4)**-0.4
        # xsh, ysh, x_rms, y_rms = flt.align_bright_objects(xspec=xspec, yspec=yspec, ds9=ds9)
        # flt.update_wcs_with_shift(xsh=xsh, ysh=ysh)
        c = flt.blot_catalog(cat, sextractor=False)
        
        ### Strip out io.fits objects
        flt.clean_for_mp()
    
    return flt
        
class GroupFLT():
    def __init__(self, files=[], refimage='Q0142-10-F140W_drz_sci.fits', segimage='Q0142-10-F140W_seg.fits', master_cat='Q0142-10-F140W.cat', pad=100):
        import multiprocessing as mp
        
        global FLTs
        
        self.FLTs = []
        self.files = files
        self.N = len(files)
        
        self.refimage = refimage
        self.segimage = segimage
        self.pad = pad
        
        self.cat = catIO.Table(master_cat)
        self.cat.rename_column('X_WORLD', 'ra')
        self.cat.rename_column('Y_WORLD', 'dec')
        self.cat.rename_column('NUMBER', 'id')
        
        #### Read in the FLTs
        self._read_FLTs(self.files, initialize=True, parallel=(mp.cpu_count() > 1))
        self.pivot = self.FLTs[0].pivot
        
        ids = [0]
        for flt in self.FLTs:
            ids = np.unique(np.append(ids, flt.catalog['id']))
        
        self.ids = ids[0:]
        
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
        
        self.dispersion_PA = np.array([flt.dispersion_PA for flt in self.FLTs])
        self.round_dispersion_PA = np.cast[int](np.round(self.dispersion_PA))
        self.grism = np.array([flt.grism for flt in self.FLTs])
        
    def _read_FLTs(self, files, initialize=True, parallel=True):
        import multiprocessing as mp
        global FLTs

        if initialize:
            FLTs = []
            
        if parallel:
            t0_pool = time.time()
            pool = mp.Pool(processes=mp.cpu_count())
    
            results = [pool.apply_async(_loadFLT, (self.files[i], self.refimage, self.segimage, self.pad, self.cat)) for i in xrange(self.N)]
            pool.close()
            pool.join()
        
            for res in results:
                flt_i = res.get(timeout=1)
                FLTs.append(flt_i)
            
            t1_pool = time.time()
        else:
            for i in xrange(self.N):
                flt_i = _loadFLT(self.files[i], self.refimage, self.segimage, self.pad, self.cat)
                FLTs.append(flt_i)
        
        ### Need to re-open the FLT files.  Don't bother with refimage/segimage 
        for flt in FLTs:
            flt.im = pyfits.open(flt.flt_file)
                
        self.FLTs = FLTs
            
    def init_fit_data(self, xspec=None, yspec=None, maglim={'A':25,'B':24,'C':22,'D':21,'E':21,'F':21}):
        if xspec is None:
            xspec = np.arange(7000,1.8e4,5)
        
        if yspec is None:
            yspec = (xspec/self.pivot)**-1 #-0.2
        
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
    
    def get_beams(self, id=995, fcontam=1, cutout_dimensions=[14,14], ds9=None, flt_mask=None, beam_name='A'):
        """
        Get beams for a given object
        """
        beams = []
        for ii, flt in enumerate(self.FLTs):
            if flt_mask is not None:
                if not flt_mask[ii]:
                    continue
                    
            if id not in flt.catalog['id']:
                continue
            
            ix = flt.catalog['id'] == id
            xi = flt.catalog['x_flt'][ix][0]
            yi = flt.catalog['y_flt'][ix][0]
            
            if (yi < flt.pad) | (yi > flt.sh_flt[0]+flt.pad) | (xi < cutout_dimensions[1]) | (xi > flt.sh_pad[1]-cutout_dimensions[1]):
                #print 'Skip:', xi, yi
                continue
            else:
                pass
                #print xi, yi
            
            beam = mywfc3.flt.model.BeamCutout(x=xi, y=yi, cutout_dimensions=cutout_dimensions, beam=beam_name, conf=flt.conf, GrismFLT=flt) #direct_flam=flt.flam, grism_flt=flt.im)
            beam.shc = beam.cutout_sci.shape
            if np.product(beam.shc) != np.product(beam.shg):
                continue
            
            beam.slx = slice(0, beam.shc[1])
            beam.cutout_seg = np.cast[np.float32](beam.get_flam_thumb(flt.seg))
            beam.id = id
            
            # import pysynphot as S
            # g = S.GaussianSource(1.e4, 1.4e4, 8)
            # beam.compute_model(beam.thumb, id=id, xspec=g.wave, yspec=g.flux)
            # beam.line_model = beam.model*1.
            
            beam.compute_model(beam.thumb, id=id, xspec=self.fit_data[id]['xspec'], yspec=self.fit_data[id]['yspec'])
            beam.contam = beam.get_cutout(flt.model) - beam.model[:,beam.slx]
            
            beam.ivar = np.cast[np.float32](1/(beam.cutout_err**2+(fcontam*beam.contam)**2))
            beam.ivar[(beam.cutout_err == 0) | (beam.cutout_dq > 0) | (beam.cutout_err > 0.1)] = 0
            
            yp, xp = np.indices(beam.shg)
            #xp = (xp.flatten()-beam.cutout_dimensions[0])*(beam.lam[1]-beam.lam[0])+beam.lam[0]
            #xp = (xp.flatten()-beam.cutout_dimensions[0])*(beam.lam[1]-beam.lam[0])
            beam.xpf = (xp[:,beam.slx].flatten()-beam.cutout_dimensions[0]*1.)/100.
            beam.xp = beam.xpf.reshape(beam.shc)
            beam.cleanf = (beam.cutout_sci-beam.contam).flatten()
            beam.clean = beam.cleanf.reshape(beam.shc)
            
            beam.bg = 0.
            
            beams.append(beam)
            
            if ds9 is not None:
                ds9.frame(ii+1)
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
                self.FLTs[i].model = self.FLTs[i].modelf.reshape(self.FLTs[i].sh_pad)
                #print FLTs[key].modelf.max()

            t1_pool = time.time()
        else:
            t0_pool = time.time()
            for i in xrange(self.N):
                out = _FLT_compute_model(i, fit_data, initialize, verbose)
            
            t1_pool = time.time()
                
        self.dt = t1_pool - t0_pool
    
    def refine_model(self, maglim=23, ds9=None, off=0.04):
        """
        Update the parent contamination model based on a simple continuum
        fit to the bright objects in the catalog
        
        offset to allow for negative backgrounds and non-negative coefficient fits
        """
        
        bright_ids = []
        mags = []
        
        for flt in self.FLTs:
            keep = flt.catalog['MAG_AUTO'] < maglim
            ids_i = flt.catalog['id'][keep]            
            mags_i = flt.catalog['MAG_AUTO'][keep]          
            for id, mag in zip(ids_i, mags_i):
                if id not in bright_ids:
                    bright_ids = np.append(bright_ids, id)
                    mags = np.append(mags, mag)
                    
        so = np.argsort(mags)
        N = len(mags)
        cutout_dimensions=[10,10]
        
        #off = 0.04
                
        for ii, id in enumerate(bright_ids[so]):
            try:
                beams = self.get_beams(id=id, cutout_dimensions=cutout_dimensions, ds9=None)
                if len(beams) == 0:
                    continue
            
                mb = mywfc3.flt.multifit.MultiBeam(beams, simple_templates=True)
            except:
                print 'id=%d Failed' %(id)
                continue
            
            mb.y += off
            
            ### zero out emission line templates
            #mb.A[:,2:] *= 0.
            #mb.line_template_hi[1]*=0
            #mb.line_template_lo[1]*=0 
            
            z = 0
            
            out, Ap, coeff = mb.fit_at_z(z=z, get_fit=True, fit_background=True)
            xfit, yfit = mb.get_1d_spec(coeff, z)
            self.update_model(id, xspec=xfit, yspec=yfit)
            
            print 'Refine: %4d/%4d %6d %.1f %s' %(ii,N, id, mags[so][ii], ' '.join(['%.1e' %(co) for co in coeff]))
            if ds9 is not None:
                ds9.view(self.FLTs[0].im_data['SCI']-self.FLTs[0].model)
    
    def save_data(self, file='multifit.npy'):
        data = {'fit_data': self.fit_data}
        backgrounds = {}
        for FLT in self.FLTs:
            if hasattr(FLT, 'fit_background_result'):
                backgrounds[FLT.flt_file] = FLT.fit_background_result
        
        data['background'] = backgrounds
        np.save(file, [data])
    
    def load_data(self, file='multifit.npy', verbose=True, recompute_model=True, apply_background=True):
        if verbose:
            print '\n==========\nLoad data: %s\n' %(file)
        
        indata = np.load(file)[0]
        if 'fit_data' in indata.keys():
            for id in indata['fit_data'].keys():
                if verbose:
                    print no_newline+'%s: %d > self.fit_data' %(file, id)
                    
                self.fit_data[id] = indata['fit_data'][id]
            
            if recompute_model:
                self.mp_compute_models(self.fit_data, initialize=True, verbose=1, parallel=True)
                
        if ('background' in indata.keys()) & apply_background:
            for flt_file in indata['background'].keys():
                pfit = indata['background'][flt_file]
                for FLT in self.FLTs:
                    if FLT.flt_file == flt_file:
                        if verbose:
                            print 'Apply background %s (p0_0=%7.4f)' %(flt_file, pfit.parameters[0])
                        
                        FLT.fit_background(pfit=pfit, apply=True)
                        
    def extract_and_fit(self, id=18629, zr=[0.68,2.4], dz=0.003, cutout_dimensions=[14,14], ds9=None, update_model=True, zsp=None):
        import copy
        
        ### Extract 
        beams = self.get_beams(id=id, cutout_dimensions=cutout_dimensions, ds9=ds9)
        mb = MultiBeam(beams)
        
        ### use trace geometry of the first beam
        ref, refwcs = mb.beams[0].make_wcs_header(mb.beams[0].model)
        refh = ref.header
        
        ## straighten trace
        #refh['CD2_1'] = 0
        
        ### drizzle all beams to trace of first beam
        header, sci, wht, ctx = mb.drizzle(beam_attr='cutout_sci', pixfrac=0.3, ref_header=refh, ds9=ds9)
        header, contam, wht, ctx = mb.drizzle(beam_attr='contam', pixfrac=0.3, ref_header=refh, ds9=None)
        if ds9 is not None:
            ds9.view(sci-contam)
        
        ## Make replacement MultiBeam list with just the drizzled beam
        driz_beam = copy.copy(mb.beams[0])
        driz_beam.clean = np.roll(sci-contam, -1, axis=1)
        driz_beam.cleanf = driz_beam.clean.flatten()
        driz_beam.thumb *= 0.
        for ib, beam in enumerate(beams):
            driz_beam.thumb += mb.beams[ib].thumb/len(beams)

        driz_beam.ivar = np.roll(wht, -1, axis=1).flatten()
        mbd = mywfc3.flt.multifit.MultiBeam([driz_beam])

        ## Fit the redshift
        zgrid, lnp, coeffs = mbd.fit_grid(zr=zr, dz=dz)
        z_gris, lnp_max = get_parabola_max(zgrid, lnp)
        
        plt.plot(zgrid, lnp-lnp_max)
        plt.scatter(z_gris, 0, color='k')

        if zsp is not None:
            print 'dz/(1+z): %.4f' %((z_gris-zsp)/(1+zsp))
            plt.plot(zsp*np.ones(2), [(lnp-lnp.max()).min(), 0], color='k')
        
        ### Get the best-fit spectrum
        out, coeff = mb.fit_at_z(z=z_gris, get_fit=True)
        xfit, yfit = mb.get_1d_spec(coeff, z_gris)

        mb.reshape_modelfit(out, attr='modelfit')
        
        ### Optimal extractions
        refh_out = refh.copy()
        refh_out['CD2_1'] = 0
        header, sci, wht, ctx = mb.drizzle(beam_attr='cutout_sci', pixfrac=0.3, ref_header=refh_out, ds9=None)
        header, contam, wht, ctx = mb.drizzle(beam_attr='contam', pixfrac=0.3, ref_header=refh_out, ds9=None)
        header, msci, mwht, mctx = mb.drizzle(beam_attr='modelfit', pixfrac=0.3, ref_header=refh_out, ds9=None)
        
        if update_model:
            self.update_model(id, xspec=xfit, yspec=yfit)
                
    def fit_multiple(self, ids=[]):
        """ Doesn't work
        """
        import multiprocessing as mp
        
        t0_pool = time.time()
        pool = mp.Pool(processes=mp.cpu_count())

        results = [pool.apply_async(func_fit, (self, id)) for id in ids]
        pool.close()
        pool.join()
    
        for res in results:
            status = res.get(timeout=1)
        
        t1_pool = time.time()
        
def func_fit(self, id):
    self.extract_and_fit(id=id)
    return 1
    
def regularize_templates():
    """
    Make sure the EAZY grism line templates all have the same length
    """
    templates = ['templates/EAZY_v1.0_lines/eazy_v1.0_sed1_nolines.dat',
    'templates/EAZY_v1.0_lines/eazy_v1.0_sed2_nolines.dat',  
    'templates/EAZY_v1.0_lines/eazy_v1.0_sed3_nolines.dat',     
    'templates/EAZY_v1.0_lines/eazy_v1.0_sed4_nolines.dat',     
    'templates/EAZY_v1.0_lines/eazy_v1.0_sed5_nolines.dat',     
    'templates/EAZY_v1.0_lines/eazy_v1.0_sed6_nolines.dat',     
    'templates/cvd12_t11_solar_Chabrier.extend.dat',     
    'templates/dobos11/SF0_0.emline.loOIII.txt',    
    'templates/dobos11/SF0_0.emline.hiOIII.txt',   
    'templates/dobos11/bc03_pr_ch_z02_ltau07.0_age09.2_av2.5.dat']
    
    for t in templates:
        tx, ty = np.loadtxt('/Users/brammer/3DHST/Spectra/Work/%s' %(t), unpack=True)
        print t, len(tx)
    
    tx_ref, ty = np.loadtxt('/Users/brammer/3DHST/Spectra/Work/%s' %(templates[-2]), unpack=True)
    for t in templates:
        tx, ty = np.loadtxt('/Users/brammer/3DHST/Spectra/Work/%s' %(t), unpack=True)
        if len(tx) != len(tx_ref):
            ty_ref = unicorn.utils_c.interp_conserve_c(tx_ref, tx, ty)
        else:
            ty_ref = ty
        
        if 'emline' in t:
            plt.plot(tx_ref, ty_ref)
            
        np.save('/Users/brammer/3DHST/Spectra/Work/templates/match_grism_templates/' + os.path.basename(t) + '.npy', [tx_ref, ty_ref/np.trapz(ty_ref, tx_ref)*1.e2])
            
class MultiBeam(object):
    def __init__(self, beams=[], model_min=0.01, ivar_min=[0.05,97], simple_templates=False):
        self.beams = beams
        import pysynphot as S
        #self.line_template_hi = np.load(os.getenv('THREEDHST')+'/templates/dobos11/SF0_0.emline.hiOIII.txt.npy')[0]        
        #self.line_template_lo = np.load(os.getenv('THREEDHST')+'/templates/dobos11/SF0_0.emline.emline.m.txt.npy')[0]
        #self.line_template_lo = np.loadtxt(os.getenv('THREEDHST')+'/templates/cvd12_t11_solar_Chabrier.extend.m.dat', unpack=True)
               
        ### SImple line templates
        if False:
            xarr = np.logspace(2,5,100)
            yarr = xarr*0
            cont_spec = S.ArraySpectrum(xarr, yarr, fluxunits='flam')
            doublet = cont_spec + 1*S.GaussianSource(1., 3729, 10) + 1*S.GaussianSource(1., 4863, 10) + 1*S.GaussianSource(1., 4959, 10) + 2.98*S.GaussianSource(1., 5007, 10) + 0*S.GaussianSource(1, 6563., 10) + 0*S.GaussianSource(1, 6725, 10) 
            self.line_template_lo = [doublet.wave, doublet.flux]
            single = cont_spec + 0*S.GaussianSource(1., 3727, 10) + 0.*S.GaussianSource(1., 4861, 10) + 0*S.GaussianSource(1., 4959, 10) + 0*S.GaussianSource(1., 5007, 10) + 1*S.GaussianSource(1, 6563., 10) + 0.2*S.GaussianSource(1, 6731, 10) 
            self.line_template_hi = [single.wave, single.flux]
        
        if simple_templates:
            if len(self.beams) > 0:
                lam = self.beams[0].lam*1
            else:
                lam = np.arange(7000,1.6e4,40)
                
            #med = np.median(lam)
            max = lam.max()*1.05
            min = lam.min()*0.95
            lam = np.arange(min, max, 10)
            width = max - min
            self.templates = np.array([(lam-min)/width, (max-lam)/width])
            self.template_wave = lam
        else:
            template_files = ['eazy_v1.0_sed1_nolines.dat',
                'eazy_v1.0_sed2_nolines.dat',  
                'eazy_v1.0_sed3_nolines.dat',     
                'eazy_v1.0_sed4_nolines.dat',     
                'eazy_v1.0_sed5_nolines.dat',     
                'eazy_v1.0_sed6_nolines.dat',     
                'cvd12_t11_solar_Chabrier.extend.dat',     
                'SF0_0.emline.loOIII.txt',    
                'SF0_0.emline.hiOIII.txt',   
                'bc03_pr_ch_z02_ltau07.0_age09.2_av2.5.dat']
        
            template_files = ['/Users/brammer/3DHST/Spectra/Work/templates/match_grism_templates/%s' %(file) for file in template_files]
            self.template_wave = np.load(template_files[0]+'.npy')[0]
            self.templates = np.array([np.load(file+'.npy')[1] for file in template_files])
            
        self.NTEMP = self.templates.shape[0]
        
        #### Initialize multifit arrays
        A = None
        y = None
        mask = None
        ivar = None
        
        for ib, beam in enumerate(beams):
            beam.compute_model(beam.thumb, id=beam.id)
            beam.cmodel = beam.model[:,beam.slx].flatten()
        
            maski = (beam.ivar.flatten() != 0) & (beam.cmodel > model_min*beam.cmodel.max())
            
            if A is None:
                A = np.zeros((beam.modelf.shape[0], self.NTEMP+len(beams)), dtype=np.float64)
                #A = np.vstack([(beam.xpf*beam.cmodel), beam.cmodel, beam.cmodel*0., beam.cmodel*0., np.ones(beam.cmodel.size)]).T
                A[:, self.NTEMP+ib] = 1.
                y = beam.cleanf*1
                mask = maski
                ivar = beam.ivar.flatten()
            else:
                Ap = np.zeros((beam.modelf.shape[0], self.NTEMP+len(beams)), dtype=np.float64)
                Ap[:,self.NTEMP+ib] = 1.
                #Ap = np.vstack([(beam.xpf*beam.cmodel), beam.cmodel, beam.cmodel*0, beam.cmodel*0., np.ones(beam.cmodel.size)]).T
                A = np.append(A, Ap, axis=0)
                y = np.append(y, beam.cleanf*1-0.0)
                mask = np.append(mask, maski)
                ivar = np.append(ivar, beam.ivar.flatten())
        
        self.A = A
        self.y = np.cast[np.float64](y)
        self.ivar = np.cast[np.float64](ivar)
        try:
            self.mask = mask & (self.ivar > ivar_min[0]*np.percentile(self.ivar[mask], ivar_min[1]))
        except:
            self.mask = mask
            
    def fit_at_z(self, z=0, coeffs=None, get_fit=False, fit_background=True, nolines=False):
        """
        TBD: add photometric constraints, at least from the direct image.
             - should help fits of bright continuum objects
             
        """
        import sklearn.linear_model
        
        Ap = self.A*1
        for itemp in range(self.NTEMP):
            for i, beam in enumerate(self.beams):
                lmodel_i = beam.compute_model(beam.thumb, id=beam.id, xspec=self.template_wave*(1+z), yspec=self.templates[itemp,:], in_place=False)
                if i == 0:
                    lmodel = lmodel_i*1
                else:
                    lmodel = np.append(lmodel, lmodel_i)
            
            Ap[:,itemp] = lmodel
                
        if coeffs is None:
            # clf = sklearn.linear_model.LinearRegression(normalize=False, fit_intercept=False)
            # status = clf.fit(Ap[self.mask,:], self.y[self.mask])
            # coeffs = clf.coef_
            import unicorn.utils_c
            if fit_background:
                ### Fit background
                mask_nneg = self.mask & (self.y > 0)
                amat = unicorn.utils_c.prepare_nmf_amatrix(self.ivar[mask_nneg], Ap[mask_nneg,:].T)
                init = np.ones(Ap.shape[1])
                init[self.NTEMP:] = 0.05
                # if hasattr(self, 'coeffs'):
                #     init=self.coeffs
                coeffs = unicorn.utils_c.run_nmf(self.y[mask_nneg]+0.05, self.ivar[mask_nneg], Ap[mask_nneg,:].T, amat, init_coeffs=init)
                coeffs[self.NTEMP:] -= 0.05
                #self.coeffs = coeffs
            else:
                ### No background
                mask_nneg = self.mask & (self.y > 0)
                amat = unicorn.utils_c.prepare_nmf_amatrix(self.ivar[mask_nneg], Ap[mask_nneg,:self.NTEMP].T)
                coeffs = unicorn.utils_c.run_nmf(self.y[mask_nneg], self.ivar[mask_nneg], Ap[mask_nneg,:self.NTEMP].T, amat, toler=1.e-5)
                coeffs = np.append(coeffs, np.zeros(len(self.beams)))
            
            if nolines:
                coeffs[-4:-2] = 0.
                
        
        out = np.dot(Ap, coeffs)
        if get_fit:
            return out, Ap, coeffs #clf.coef_
            
        lnp = -0.5*np.sum((out-self.y)[self.mask]**2*self.ivar[self.mask])
        
        return lnp, coeffs #clf.coef_
    
    def get_1d_spec(self, coeffs, z):
        #beam = self.beams[0]
        #xpix = (self.line_template_hi[0]*(1+z) - beam.lam[0])/(beam.lam[1]-beam.lam[0])/100.
        #yspec = (coeff[0]*xpix+coeff[1]) + self.line_template_hi[1]*coeff[2] + self.line_template_lo[1]*coeff[3]
        #return self.line_template_hi[0]*(1+z), yspec
        return self.template_wave*(1+z), np.dot(coeffs[:self.NTEMP], self.templates) #+ coeffs[-1]
        
    def reshape_modelfit(self, out, attr='modelfit'):
        i0 = 0
        for beam in self.beams:
            i1 = beam.modelf.size
            beam.__setattr__(attr, out[i0:i0+i1].reshape(beam.shg))
            i0 += i1
    
    def get_zgrid(self, zr=[0.7,3.4], dz=0.01):
        zgrid = np.exp(np.arange(np.log(1+zr[0]), np.log(1+zr[1]), dz))-1
        return zgrid
        
    def fit_grid(self, zr=[0.7, 3.4], dz=0.01, zgrid = None, verbose=True, fit_background=True):
        """
        """
                
        if zgrid is None:
            zgrid = self.get_zgrid(zr, dz) # np.exp(np.arange(np.log(1+zr[0]), np.log(1+zr[1]), dz))-1
        
        nz = len(zgrid)
        lnp = np.zeros(nz)
        coeffs = np.zeros((nz, self.A.shape[1]))
        for i in range(nz):
            lnp[i], coeffs[i,:] = self.fit_at_z(z=zgrid[i], fit_background=fit_background)
            if verbose:
                print no_newline + 'id=%d: %.4f  - %.4f' %(self.beams[0].id, zgrid[i], lnp[i])
                
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
        
    def drizzle(self, beam_attr='clean', ref_header=None, pixfrac=1, ds9=None, update_fit=True, kernel='square', thumb=False):
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
            if thumb:
                beam_wcs.wcs.crval[0] = ref_wcs.wcs.crval[0]
                beam_wcs.wcs.cd[0,:] = ref_wcs.wcs.cd[0,:]
                
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
        if thumb:
            return outsci[:, :outsci.shape[0]]
        else:
            return ref_header, outsci, outwht, outctx
        
        ########### DOne
        
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

def manybeams():
    import glob
    import multiprocessing as mp
    import time
    
    files=glob.glob('*[0-9].npy')[:14]

    t0 = time.time()
    
    pool = mp.Pool(processes=4)
    results = [pool.apply_async(fit_MultiBeam, (file,)) for file in files]
    pool.close()
    #pool.join()
    
    while not results[-1]._ready:
        print time.time()
        for i, res in enumerate(results):
            print i, res._ready
        pass
    
    print 'Done?'
        
    status = res.get(timeout=1)
        
    t1 = time.time()
    
    t0x = time.time()
    for file in files:
        fit_MultiBeam(file)

    t1x = time.time()
      
    print 'Parallel: %.2f' %(t1-t0)  
    print 'Serial  : %.2f' %(t1x-t0x)  

def fit_MultiBeam(file='udf_multifit_21798.npy'):
    import copy
    import scipy.ndimage as nd
    
    print '%s start' %(file)
    
    mb = np.load(file)[0]
    
    ref, refwcs = mb.beams[0].make_wcs_header(mb.beams[0].model)
    refh = ref.header

    header, sci, wht, ctx = mb.drizzle(beam_attr='cutout_sci', pixfrac=0.3, ref_header=refh, ds9=None)
    header, contam, wht, ctx = mb.drizzle(beam_attr='contam', pixfrac=0.3, ref_header=refh, ds9=None)
    #ds9.view(sci-contam)

    driz_beam = copy.copy(mb.beams[0])
    driz_beam.clean = np.roll(sci-contam, -1, axis=1)
    driz_beam.cleanf = driz_beam.clean.flatten()
    driz_beam.thumb *= 0.
    for ib, beam in enumerate(mb.beams):
        driz_beam.thumb += mb.beams[ib].thumb/len(mb.beams)
    
    #driz_beam.cleanf = (sci-contam).flatten()
    driz_beam.ivar = np.roll(wht, -1, axis=1).flatten()
    
    #### Fit drizzled beam
    mbd = mywfc3.flt.multifit.MultiBeam([driz_beam])
    
    zr = [0.7,2.4]
    dz = 0.002
    
    try:
        zgrid, lnp, coeffs = mbd.fit_grid(zr=zr, dz=dz)
    except:
        return False
        
    #lnp = np.log(nd.gaussian_filter(np.exp(lnp-lnp.max()), 2))
    lnp = nd.gaussian_filter(lnp, dz/0.001)
    
    try:
        z_gris, lnp_max = get_parabola_max(zgrid, lnp)
    except:
        return False
        
    out, Ap, coeff = mb.fit_at_z(z=z_gris, get_fit=True)
    xfit, yfit = mb.get_1d_spec(coeff, z_gris)
    #group.update_model(id, xspec=xfit, yspec=yfit)
    
    mb.reshape_modelfit(out, attr='modelfit')
    #refh['CRPIX2'] += 0.1
    header, msci, mwht, mctx = mb.drizzle(beam_attr='modelfit', pixfrac=0.3, ref_header=refh)
    
    wave = (np.arange(refh['NAXIS1'])-refh['CRPIX1'])*refh['CD1_1'] + refh['CRVAL1']
    
    print '%s done, z_gris=%.3f' %(file, z_gris)
    
    np.save(file.replace('.npy','.fit.npy'), [driz_beam, msci, wave, z_gris, xfit, yfit, out, coeff])
    return z_gris, xfit, yfit
    
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
                    #if (beam == 'F') & (flt.grism != 'G102'):
                    if beam not in flt.conf.beams:
                        continue
                    
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
    
def _test_figure():
    groups = [g102, g141]
    driz_list = []
    for g in groups:
        full_beams = g.get_beams(id=id, cutout_dimensions=cutout_dimensions, ds9=None)
        for beams in [full_beams]:#, full_beams[8:]]:
            mb = mywfc3.flt.multifit.MultiBeam(beams)
    
            #refh, refwcs = mb.get_drizzle_wcs(center_wave=14000.0, dlam=46, NX=100, spatial_scale=1, NY=cutout_dimensions[0])
            #refh['CRPIX2'] -= 2

            ref, refwcs = mb.beams[0].make_wcs_header(mb.beams[0].model)
            refh = ref.header
            #refh['CRVAL2'] *= 2# 1
            #refh['CD2_1'] = 0
    
            header, sci, wht, ctx = mb.drizzle(beam_attr='cutout_sci', pixfrac=0.3, ref_header=refh, ds9=ds9b)
            #header, sci, wht, ctx = mb.drizzle(beam_attr='clean', pixfrac=0.3, ref_header=refh, ds9=ds9b)
            header, contam, cwht, ctx = mb.drizzle(beam_attr='contam', pixfrac=0.3, ref_header=refh, ds9=None)
            header, fmsci, mwht, mctx = mb.drizzle(beam_attr='model', pixfrac=0.3, ref_header=refh)
            ds9b.view(sci-contam)
    
            driz_beam = copy.copy(mb.beams[0])
            driz_beam.cutout_sci = np.roll(sci, -1, axis=1)
            driz_beam.contam = np.roll(contam, -1, axis=1)
            driz_beam.clean = np.roll(sci-contam, -1, axis=1)
            driz_beam.cleanf = driz_beam.clean.flatten()
            driz_beam.thumb *= 0.
            for ib, beam in enumerate(beams):
                driz_beam.thumb += mb.beams[ib].thumb/len(beams)
        
            #driz_beam.cleanf = (sci-contam).flatten()
            min_err = 0.05
            wht = 1./(1/wht+(min_err*fmsci)**2)
            driz_beam.ivar = np.roll(wht, -1, axis=1).flatten()
            driz_list.append(driz_beam)
        
    mbd = mywfc3.flt.multifit.MultiBeam(driz_list)
    #mbd = mb
    
    zgrid, lnp, coeffs = mbd.fit_grid(zr=zr, dz=dz)
    #lnp = np.log(nd.gaussian_filter(np.exp(lnp-lnp.max()), 2))
    lnp = nd.gaussian_filter(lnp, dz/0.001)
    
    z_gris, lnp_max = get_parabola_max(zgrid, lnp)
    #plt.plot(zgrid, lnp-lnp_max)
    #plt.scatter(z_gris, 0, color='k')
    
    if zsp:
        print 'dz/(1+z): %.4f' %((z_gris-zsp)/(1+zsp))
        #plt.plot(zsp*np.ones(2), [(lnp-lnp.max()).min(), 0], color='k')
    
    out, Ap, coeff = mbd.fit_at_z(z=z_gris, get_fit=True)
    xfit, yfit = mbd.get_1d_spec(coeff, z_gris)
    #group.update_model(id, xspec=xfit, yspec=yfit)
    
    mbd.reshape_modelfit(out, attr='modelfit')
    #refh['CRPIX2'] += 0.1
    #header, msci, mwht, mctx = mb.drizzle(beam_attr='modelfit', pixfrac=0.3, ref_header=refh)
    
    import matplotlib.gridspec as gridspec
    
    gs = gridspec.GridSpec(6, len(mbd.beams))
    gst = gridspec.GridSpec(6, len(mbd.beams))
    gst.update(hspace=0.1)
    
    fig = plt.figure(figsize=[5*len(mbd.beams),7])
    
    spec_axes, all_max = [], 0
    
    for ib, beam in enumerate(mbd.beams):
        ref, refwcs = beam.make_wcs_header(beam.model)
        refh_out = ref.header
        wave = (np.arange(refh_out['NAXIS1'])-refh_out['CRPIX1'])*refh_out['CD1_1'] + refh_out['CRVAL1']
    
        ### Optimal extractions
        beam.compute_model(beam.thumb, id=beam.id)
        prof = beam.model/beam.model.sum(axis=0)
        
        owht = beam.ivar.reshape(beam.shg)
        num = prof*((beam.cutout_sci-coeff[mbd.NTEMP+ib])*owht)
        cnum = prof*(beam.contam*owht)
        den = prof**2*owht
        opt = num.sum(axis=0)/den.sum(axis=0)
        copt = cnum.sum(axis=0)/den.sum(axis=0)
        opt_var = 1./den.sum(axis=0)
        mnum = prof*((beam.modelfit-coeff[mbd.NTEMP+ib])*owht)
        mopt = mnum.sum(axis=0)/den.sum(axis=0)
        
        ys = np.interp(wave, beam.lam, beam.ysens)
        ys[ys < 0.1*ys.max()] = np.nan
        for arr in [opt, copt, mopt]:
            arr *= 100/ys
        
        opt_var *= (100/ys)**2
        
        #max = np.percentile(osci-ocontam, 100)
        max = (beam.clean[beam.model > 0]).max()
        if wave.max() > 1.4e4:
            xr = [1.02e4,1.76e4]
            xi = np.arange(1.1,1.71,0.1)
        else:
            xr = [0.76e4,1.19e4]
            xi = np.arange(0.8,1.11,0.1)
        
        for ii, data in enumerate([beam.clean, beam.modelfit]):
            ax = fig.add_subplot(gst[ii,ib])
            ax.imshow(data, vmin=-0.1*max, vmax=max, interpolation='Nearest', aspect='auto', origin='lower')
            ax.set_xlim(np.interp(xr, wave, np.arange(beam.shg[1])))
            ax.set_ylim(beam.shg[0]/2.+1+np.array([-1,1])*beam.cutout_dimensions[0])
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.set_xticks(np.interp(xi*1.e4, wave, np.arange(beam.shg[1])))
        
        ax = fig.add_subplot(gs[2:4,ib])
        
        ax.plot(wave/1.e4, opt, alpha=0.5, color='0.5', linestyle='steps-mid')
        ax.plot(wave/1.e4, opt-copt, alpha=0.5, color='k', linestyle='steps-mid')
        ax.fill_between(wave/1.e4, opt-copt, wave*0., alpha=0.1, color='0.5')
        ax.plot(wave/1.e4, mopt, alpha=0.5, color='r', linestyle='steps-mid', linewidth=2)
        ax.errorbar(wave/1.e4, opt-copt, np.sqrt(opt_var), linestyle='None', marker='.', color='k', ecolor='k', alpha=0.5)
        ax.set_xlim(np.array(xr)/1.e4)
        #ax.set_xticklabels([]); #ax.set_yticklabels([])
        ax.set_ylabel('flux (e/s)')
        ax.set_xlabel(r'$\lambda$')
    
        if wave.max() > 1.4e4:
            ymax = (opt-copt+np.sqrt(opt_var))[(wave > 1.1e4) & (wave < 1.65e4) & np.isfinite(opt_var)].max()
        else:
            ymax = (opt-copt+np.sqrt(opt_var))[(wave > 0.85e4) & (wave < 1.12e4) & np.isfinite(opt_var)].max()
        
        ax.set_ylim(-0.3*ymax, 1.1*ymax)
        all_max = np.maximum(ymax, all_max)
        spec_axes.append(ax)
    
    for axx in spec_axes:
        axx.set_ylim(-0.15*all_max, 1.1*all_max)
        axx.grid(alpha=0.2)
    
    axx.set_yticklabels([])
        
    #### p(z)
    #gs = gridspec.GridSpec(6, 1)
    ax = fig.add_subplot(gs[4:,0])
    ax.plot(zgrid, lnp-lnp_max, color='b')
    ax.fill_between(zgrid, lnp-lnp_max, (lnp-lnp_max).min(), color='b', alpha=0.1)
    
    if zsp:
        #ax.vlines(np.array([3727, 4959, 5007, 6563])*(1+zsp)/1.e4, (lnp-lnp_max).min(), 0, color='r', alpha=0.8)
        ax.vlines([zsp], (lnp-lnp_max).min(), 0, color='r', alpha=0.8)
    
    ### 
    # eazypz = eazy.getEazyPz(np.arange(len(c))[c['id'] == id][0], MAIN_OUTPUT_FILE=os.path.basename(z.filename.split('.zout')[0]), OUTPUT_DIRECTORY=os.path.dirname(z.filename), CACHE_FILE='Same', binaries=None, get_prior=False)
    # eazy_lnp = np.log(eazypz[1])
    # ax.plot(eazypz[0], eazy_lnp - eazy_lnp.max(), color='k', alpha=0.4)
    # 
    ax.scatter(z_gris, 0)
    ymin = (lnp-lnp_max).min()
    ax.text(z_gris, -0.08*ymin, '%.3f' %(z_gris), ha='center', va='top', fontsize=9)
    ax.set_xlim(zgrid.min(), zgrid.max()) 
    #ax.set_xlim(np.array(xr)/1.e4)
   
    ax.set_xlabel(r'z'); ax.set_ylabel(r'ln $p = -0.5\chi^2$')
    #pzline = np.log(np.maximum(coeffs[:,3], 1.e-20))+lnp
    #ax.plot((1+zgrid)*5007./1.e4, pzline-pzline.max(), color='green')

    ax.set_ylim(ymin, -0.1*ymin)
        
    if len(mbd.beams) > 1:
        ax = fig.add_subplot(gs[4:,1])    
        pz = np.exp(lnp-lnp_max); pt = np.trapz(pz, zgrid)
        ax.plot(zgrid, pz, color='b')
        ax.fill_between(zgrid, pz, color='b', alpha=0.1)
        ax.set_xlabel(r'z'); ax.set_ylabel(r'$p(z)$')
        ax.set_yticklabels([])
        
    fig.tight_layout(pad=0.5)
    