
import glob
import stwcs
import numpy as np
from scipy.interpolate.fitpack import splev, splrep

import threedhst.dq
ds9 = threedhst.dq.myDS9()

import mywfc3.flt.model
import astropy.io.fits as pyfits
from threedhst import catIO

files=glob.glob('ico205l[orux]q*flt.fits')
refimage='wfi2026-4536-05-235-f140w_drz_sci.fits'

### Refsdal
info = catIO.Table('info')
pa = np.cast[int](info['PA_V3'])
pas = np.unique(pa[info['FILTER'] == 'G141'])
files=info['FILE'][(info['FILTER'] == 'G141') & (pa == 111)]
files=info['FILE'][(info['FILTER'] == 'G141') & (pa == 119)]
#files=info['FILE'][(info['FILTER'] == 'G141') & (np.abs(pa-115) < 5)]
#files=info['FILE'][(info['FILTER'] == 'G141')]
refimage = '../MACS1149/Catalog/MACS1149-F140W_drz_sci.fits'

#files=files[:1]

FLT = {}
for i, file in enumerate(files):
    print '%d/%d %s' %(i+1, len(files), file)
    g = mywfc3.flt.model.GrismFLT(file=file, refimage=refimage)
    FLT[file] = g

### Full model
import pysynphot as S
xarr = np.arange(8000,1.9e4,40)
beta = -0.5
yarr = (xarr/1.4e4)**beta*10
bp = S.ObsBandpass('wfc3,ir,f140w')
spec = S.ArraySpectrum(xarr, yarr, fluxunits='flam').renorm(1, 'flam', bp)

for i, file in enumerate(files):
    #status = self.compute_model(x=281.13, y=856.8, sh=[20,20])
    print '%d/%d %s\n' %(i+1, len(files), file)
    self = FLT[file]
    self.full_model*=0
    skip = 40
    for xs in range(skip, 901, skip):
        for ys in range(skip, 976, skip):
            print '\x1b[1A\x1b[1M' + '%d %d' %(xs, ys)
            for beam in 'ABCD':
                status = self.compute_model(x=xs, y=ys, sh=[skip/2,skip/2], beam=beam, xspec=spec.wave, yspec=spec.flux)
    
    self.m = self.full_model.reshape((1014,1014))*1.
    ds9.view(self.im['SCI'].data-self.m)
    
    #status = self.compute_model(x=507, y=507, sh=[100,100])
    
ra, dec =  306.5641, -45.61369
dim = [20,20]
beta = -1.5

ra, dec =  306.5393, -45.608193 # star
#ra, dec =  306.5547, -45.608779 # galaxy
ra, dec =  306.53132, -45.596223 # with line, z=1.140
ra, dec = 306.5249, -45.614889 # with line, z~1.0
ra, dec = 306.56274, -45.612287 # late star
ra, dec = 306.55402, -45.594289  # galaxy

### MACS1149
ra, dec = 177.40697, 22.407434 # line emitter
ra, dec = 177.38996, 22.412696 # JD
ra, dec = 177.38236, 22.405779 # line

#dim=[100,100]
dim=[50,50]
#dim = [20,20]

beams = {}
x0, y0 = None, None
for file in files:
    g = FLT[file]
    x, y = g.all_world2pix(ra, dec)
    if x0 is None:
        x0, y0 = x, y
        
    print '%s %.2f %.2f' %(file, x-x0, y-y0)
    
    beam = mywfc3.flt.model.BeamCutout(x=x, y=y, conf=g.conf, cutout_dimensions=dim)
    gris = pyfits.open(file)
    beam.cutout = beam.get_cutout(gris['SCI'].data)
    beam.cutout_m = beam.get_cutout(g.m)
    beam.cutout_dq = beam.get_cutout(gris['DQ'].data)
    beam.hdu = beam.make_wcs_header(data=beam.cutout)
    beam.thumb = beam.get_flam_thumb(g.flam)
    #beta = -1.5
    beam.compute_model(beam.thumb, xspec=spec.wave, yspec=spec.flux)
    beam.omodel = beam.model*1.
    
    ds9.view((beam.cutout-beam.cutout_m+beam.model)*(beam.cutout_dq == 0))
    beams[file] = beam
    
#
#beta = -1.5
beam.compute_model(beam.thumb, xspec=beam.lam, yspec=(beam.lam/1.4e4)**beta)
mspec = beam.modelf.reshape(beam.shg).sum(axis=0)
#spec = beam.cutout.sum(axis=0)

### Try drizzling
NY, NX = 60, 360
data = np.zeros((NY, NX))
hdu = pyfits.HDUList([pyfits.ImageHDU(data=data)])

h = hdu[0].header

h['CRPIX1'] = NX/2
h['CRVAL1'] = 1.4e4
h['CD1_1'] = 22.5
h['CD1_2'] = 0

h['CRPIX2'] = NY/2
h['CRVAL2'] = np.interp(1.4e4, beam.lam, beam.dy)
h['CD2_1'] = 0.0#1
h['CD2_2'] = 0.5

h['CD1_1'] = 46
h['CD2_2'] = 1

#h['CD1_1'] = 92
#h['CD2_2'] = 2



h2 = h.copy()
h2['CUNIT1'] = 'Angstrom'
h2['CTYPE1'] = 'WAVE'

out_wcs = stwcs.wcsutil.HSTWCS(hdu, ext=0)

sci = np.zeros((NY, NX), dtype=np.float32)
wht = np.zeros((NY, NX), dtype=np.float32)
msci = np.zeros((NY, NX), dtype=np.float32)
mwht = np.zeros((NY, NX), dtype=np.float32)
csci = np.zeros((NY, NX), dtype=np.float32)
cwht = np.zeros((NY, NX), dtype=np.float32)
ctx = np.zeros((NY, NX), dtype=np.int16)

lbsci = np.zeros((NY, NX), dtype=np.float32)
lbwht = np.zeros((NY, NX), dtype=np.float32)

for file in files:
    # sci = np.zeros(beam.cutout.shape, dtype=np.float32)
    # wht = np.zeros(beam.cutout.shape, dtype=np.float32)
    # ctx = np.zeros(beam.cutout.shape, dtype=np.int16)

    beam = beams[file]
    
    mask = beam.cutout_dq == 0
    in_wht = np.cast[np.float32](mask)
    
    sh = beam.model.shape
    #c1 = beam.cutout*1; c1[:,sh[0]/2:-sh[0]/2] /= beam.ysens
    #c2 = beam.model*1; c2[:,sh[0]/2:-sh[0]/2] /= beam.ysens
    c1 = beam.cutout*1
    c2 = beam.model*1
    c3 = beam.cutout_m*1
    astrodrizzle.adrizzle.do_driz(c1, beam.hdu[1], in_wht, out_wcs, sci, wht, ctx, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=0.6, kernel='square', fillval=0, stepsize=10, wcsmap=None)   
    astrodrizzle.adrizzle.do_driz(c2, beam.hdu[1], in_wht, out_wcs, msci, mwht, ctx, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=0.6, kernel='square', fillval=0, stepsize=10, wcsmap=None)   
    astrodrizzle.adrizzle.do_driz(c3, beam.hdu[1], in_wht, out_wcs, csci, cwht, ctx, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=0.6, kernel='square', fillval=0, stepsize=10, wcsmap=None)   
    
    llmodel = beam.compute_model(beam.clip_thumb, xspec=rn.wave, yspec=rn.flux, in_place=False).reshape(beam.shg)
    astrodrizzle.adrizzle.do_driz(llmodel, beam.hdu[1], in_wht, out_wcs, lbsci, lbwht, ctx, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=0.6, kernel='square', fillval=0, stepsize=10, wcsmap=None)   
    
    ds9.frame(1)
    #ds9.view(sci/1.e9, header=h2)
    #ds9.frame(2)
    ds9.view((sci-csci)/1.e9, header=h2)
    
wave = (np.arange(h['NAXIS1'])-h['CRPIX1'])*h['CD1_1']+h['CRVAL1']
sci1 = sci[28:33,:].sum(axis=0)
msci1 = msci[28:33,:].sum(axis=0)
csci1 = csci[28:33,:].sum(axis=0)
lbsci1 = lbsci[28:33,:].sum(axis=0)

    #time.sleep(2)

mask = msci > 0.01*msci[np.isfinite(msci)].max()
plt.plot(((sci-csci+msci)*mask).sum(axis=0), color='black', linestyle='steps-mid')
plt.plot((msci*mask*1.).sum(axis=0), color='red', linestyle='steps-mid', linewidth=2)
plt.plot(((csci-msci)*mask*1.).sum(axis=0), color='orange', linestyle='steps-mid', linewidth=2)


### Def try fitting
def setup_fit():
    
    snlim = 0
    snlim = 1./np.sqrt(len(beams))
    
    Rmax = 5
    
    for key in FLT.keys():
        beam = beams[key]
        beam.cutout_sci = beam.get_cutout(FLT[key].im['SCI'].data)*1
        beam.cutout_dq = beam.get_cutout(FLT[key].im['DQ'].data)*1
        beam.cutout_err = beam.get_cutout(FLT[key].im['ERR'].data)*1
        
        beam.cutout_mask = (beam.cutout_dq-(beam.cutout_dq & 512) == 0) & (beam.cutout_err > 0) & (beam.cutout_err < 10)
        
        sh = beam.cutout_sci.shape
        yp, xp = np.indices(beam.thumb.shape)
        r = np.sqrt((xp-sh[0]/2)**2+(yp-sh[0]/2)**2)
        beam.norm = beam.thumb[r <= Rmax].sum()/1.e-17
        beam.clip_thumb = beam.thumb*(r <= Rmax)
        beam.compute_model(beam.clip_thumb)
        
        sn = beam.model/beam.cutout_err
        #beam.mask &= (sn > 0.5) & (beam.err > 0)
        beam.cutout_mask &= (sn > snlim) & (beam.cutout_err > 0)
        
        beam.cutout_sci[~beam.cutout_mask] = 0
        beam.cutout_err[~beam.cutout_mask] = 0
        beam.cutout_var = beam.cutout_err**2
        
        beam.cutout_scif = beam.cutout_sci.flatten()
        beam.cutout_varf = beam.cutout_var.flatten()
        beam.cutout_maskf = beam.cutout_mask.flatten()
                
        ds9.view((beam.cutout_scif-beam.cutout_modelf).reshape(beam.cutout_sci.shape)*beam.cutout_mask)
        
    
    ##
    x_nodes = np.arange(1.0e4,1.71e4,0.1e4)
    x_nodes = np.arange(1.0e4,1.71e4,0.03e4)
    x_nodes = np.arange(1.05e4,1.71e4,0.075e4)
    
    init = np.append([0, 1.148, 500], np.ones(len(x_nodes)))
    init = np.append([0, 1.2078, 500], np.ones(len(x_nodes)))
    step_sig = np.append([0.01,0.1,50], np.ones(len(x_nodes))*0.1)    

    init = np.append([0, 1.0, 500], np.ones(len(x_nodes)))
    step_sig = np.append([0.0,0.2,150], np.ones(len(x_nodes))*0.1)    

    ### don't fit redshift
    #init = np.append([0, 1.0, 0], np.ones(len(x_nodes)))
    #step_sig = np.append([0.0,0.,0], np.ones(len(x_nodes))*0.2)    
    
    obj_fun = _obj_beam_fit
    obj_args = [beams, x_nodes]
        
    ndim, nwalkers = len(init), len(init)*2
    p0 = [(init+np.random.normal(size=ndim)*step_sig) 
          for i in xrange(nwalkers)]
    
    NTHREADS, NSTEP = 2, 5
    sampler = emcee.EnsembleSampler(nwalkers, ndim, obj_fun, args = obj_args, 
                                    threads=NTHREADS)
    #
    t0 = time.time()
    result = sampler.run_mcmc(p0, NSTEP)
    t1 = time.time()
    print 'Sampler: %.1f s' %(t1-t0)
    
    param_names = ['bg', 'z', 'Ha']
    param_names.extend(['spl%d' %i for i in range(len(init)-3)])
    
    import unicorn.interlace_fit
    chain = unicorn.interlace_fit.emceeChain(chain=sampler.chain, param_names=param_names)
    
    #obj_fun(chain.map, beams)
    #obj_fun(init, beams, x_nodes)
    params = chain.map*1.
    #params = init
    
    ### Show spectra
    from scipy.interpolate.fitpack import splev, splrep

    # NDRAW=100
    # draw = chain.draw_random(NDRAW)
    # for i in range(NDRAW):
    #     line = S.GaussianSource(draw[i,2], 6563.*(1+draw[i,1]), 10)
    #     tck = splrep(x_nodes, draw[i,3:], k=3, s=0)
    #     xcon = np.arange(0.9e4,1.8e4,0.01e4)
    #     ycon = splev(xcon, tck, der=0, ext=0)
    #     spec = S.ArraySpectrum(xcon, ycon, fluxunits='flam', keepneg=True)+line
    #     plt.plot(spec.wave, spec.flux, alpha=0.1, color='red')
    #     
    line = S.GaussianSource(params[2], 6563.*(1+params[1]), 10)
    tck = splrep(x_nodes, params[3:], k=3, s=0)
    xcon = np.arange(0.9e4,1.8e4,0.01e4)
    ycon = splev(xcon, tck, der=0, ext=0)
    pspec = S.ArraySpectrum(xcon, ycon, fluxunits='flam', keepneg=True)+line

    for key in beams.keys():
        beam = beams[key]
        beam.compute_model(beam.clip_thumb, xspec=pspec.wave, yspec=pspec.flux, in_place=True)
    
    xfull, yfull, mfull, cfull = [], [], [], []
    for i, key in enumerate(beams.keys()):
        beam = beams[key]
        #ds9.frame(i+1)
        #ds9.view((beam.sci-beam.model)*beam.mask)
                
        ls = 'steps-mid'
        marker = 'None'
        #ls = '-'
        #ls, marker = 'None', '.'
        xfull = np.append(xfull, beam.lam)
        ybeam = (beam.cutout_sci-params[0]*beam.cutout_mask).sum(axis=0)[sh[0]/2:-sh[0]/2]/beam.ysens/beam.norm
        yfull = np.append(yfull, ybeam)
        #plt.plot(beam.lam, ybeam, color='black', alpha=0.5, linestyle=ls, marker=marker)
        cbeam = ((beam.cutout_m-beam.omodel)*beam.cutout_mask).sum(axis=0)[sh[0]/2:-sh[0]/2]/beam.ysens/beam.norm
        cfull = np.append(cfull, cbeam)
        mbeam = ((beam.model)*beam.mask).sum(axis=0)[sh[0]/2:-sh[0]/2]/beam.ysens/beam.norm
        mfull = np.append(mfull, mbeam)        
        plt.plot(beam.lam, ybeam-cbeam, color='black', linestyle=ls, alpha=0.5)
        plt.plot(beam.lam, mbeam, color='blue', linestyle=ls, alpha=0.5)
        
    #nx = 20
    kern = np.ones(nx)/nx
    so = np.argsort(xfull)
    plt.plot(xfull[so][nx/2::nx], nd.convolve((yfull-cfull)[so], kern)[nx/2::nx], color='orange', linewidth=2, alpha=0.8, linestyle='steps-mid')
    plt.plot(xfull[so][nx/2::nx], nd.convolve(mfull[so], kern)[nx/2::nx], color='green', linewidth=2, alpha=0.8, linestyle='steps-mid')
    
    plt.ylim(0,5)
    plt.xlim(1.e4,1.78e4)
    
    from scipy.optimize import fmin_bfgs
    p0 = fmin_bfgs(_loss, init[3:], args=(beams, x_nodes), gtol=1.e-3, epsilon=1.e-3, maxiter=100)
    params = np.array(init)
    params[3:] = p0
    
def _obj_beam_fit(params, beams, x_nodes):
    """
    params = [2.953e-03, 1.156e+00, 1.297e+01, 9.747e-01 ,  9.970e-01 ,  8.509e-01, 1.076e+00 ,  1.487e+00 ,  7.864e-01,   1.072e+00 ,  1.015e+00]
    """
    import time
    from scipy.interpolate.fitpack import splev, splrep
    import pysynphot as S
    
    ### Spline continuum + gaussian line
    l0 = 6563.*(1+params[1])
    if (l0 < 1.12e4) | (l0 > 1.63e4):
        return -np.inf
        
    line = S.GaussianSource(params[2], l0, 10)
    tck = splrep(x_nodes, params[3:], k=3, s=0)
    xcon = np.arange(0.9e4,1.8e4,0.01e4)
    ycon = splev(xcon, tck, der=0, ext=0)
    spec = S.ArraySpectrum(xcon, ycon, fluxunits='flam', keepneg=True)+line
    
    lnprob = 0
    for key in beams.keys():
        beam = beams[key]
        modelf = beam.compute_model(beam.clip_thumb, xspec=spec.wave, yspec=spec.flux, in_place=False)
        lnprob += -0.5*np.sum(((beam.cutout_scif-params[0]-modelf)**2/beam.cutout_varf)[beam.cutout_maskf])
    
    if ~np.isfinite(lnprob):
        lnprob = -np.inf
        
    print params, lnprob
    #time.sleep(0.2)
    
    return lnprob
    
def _loss(params, beams, x_nodes):
    import time
    from scipy.interpolate.fitpack import splev, splrep
    import pysynphot as S
    
    ### Spline continuum + gaussian line
    # l0 = 6563.*(1+params[1])
    # if (l0 < 1.12e4) | (l0 > 1.63e4):
    #     l0 = 1.e4
    #     
    # line = S.GaussianSource(params[2], l0, 10)
    tck = splrep(x_nodes, params, k=3, s=0)
    xcon = np.arange(0.9e4,1.8e4,0.01e4)
    ycon = splev(xcon, tck, der=0, ext=0)
    spec = S.ArraySpectrum(xcon, ycon, fluxunits='flam', keepneg=True)#+line
    
    lnprob = 0
    dof = 0
    for key in beams.keys():
        beam = beams[key]
        modelf = beam.compute_model(beam.clip_thumb, xspec=spec.wave, yspec=spec.flux, in_place=False)
        lnprob += np.sum(((beam.cutout_scif-modelf)**2/beam.cutout_varf)[beam.cutout_maskf])
        dof += beam.cutout_maskf.sum()
        
    print params, lnprob, lnprob/(dof-len(params))
    #time.sleep(0.2)
    
    return lnprob
 
### Try multiprocessing
def compute_models(beams):
    """
    more beams
    
    for key in beams.keys():
        beams[key+'x'] = beams[key]
        
    """
    for key in beams:
        beam = beams[key]
        beam.compute_model(beam.thumb)

def _go_compute_model(beam):
    beam.compute_model(beam.thumb)
            
def mp_compute_models(beams):
    import multiprocessing as mp
    
    p = mp.Pool()
    results=[p.apply_async(_go_compute_model, beams[key]) for key in beams]
    p.close()
    
        
