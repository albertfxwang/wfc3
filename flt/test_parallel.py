"""
Parallel processing for FLT modeling
"""

import glob
import os
import numpy as np
import scipy.ndimage as nd
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

import drizzlepac
import drizzlepac.updatehdr

import threedhst
import threedhst.prep_flt_astrodrizzle as prep
from threedhst import catIO
import threedhst.dq

import unicorn

#### Global
FLTs = None

def mp_compute_models():
    """
    Testing multiprocessing on multiple FLTs
    """
    import multiprocessing as mp
    import multiprocessing.queues

    import mywfc3.flt
    import mywfc3.reprocess_wfc3
    
    import time
    
    global FLTs
    FLTs = {}
    
    keys = range(8)
    
    ### testing pickling
    if False:
        self = mywfc3.flt.model.GrismFLT(file='icka01t7q_flt.fits', refimage='F160W_mosaic.fits', segimage='F160W_seg.fits')
        csex = catIO.Table('/Users/brammer/3DHST/Spectra/Work/3DHST_Detection/GOODS-N_IR.cat')
        cat = self.blot_catalog(csex, sextractor=True)
        del(self.refimage)
        del(self.segimage)
        del(self.im)
        np.save('flt.npy', [self])
        self2 = np.load('flt.npy')[0]
    
    for key in keys:
        if os.path.exists('flt.npy'):
            FLTs[key] = np.load('flt.npy')[0]
        else:
            self = mywfc3.flt.model.GrismFLT(file='icka01t7q_flt.fits', refimage='F160W_mosaic.fits', segimage='F160W_seg.fits')
            #xspec = np.arange(1.e4,1.8e4)
            #yspec = (xspec/1.6e4)**-0.4
            #xsh, ysh, x_rms, y_rms = self.align_bright_objects(xspec=xspec, yspec=yspec, ds9=ds9)
            #self.update_wcs_with_shift(xsh=xsh, ysh=ysh)
            csex = catIO.Table('/Users/brammer/3DHST/Spectra/Work/3DHST_Detection/GOODS-N_IR.cat')
            cat = self.blot_catalog(csex, sextractor=True)
            FLTs[key] = self
        
    ### For testing
    for key in keys:
        FLTs[key].flam *= (key+1)
        FLTs[key].clip *= (key+1)
    
    
    ######## Queued process, has awkward I/O for getting data out
    t0_queue = time.time()
    
    queue = mp.Queue()    
    processes = [mp.Process(target=_go_compute_model, args=(key,), kwargs={'queue':queue}, name='FLT%d' %(key)) for key in keys]

    # Run processes
    for p in processes:
        p.start()

    # Exit the completed processes
    for p in processes:
        p.join()
    
    # Get the outputs
    results = [queue.get(timeout=1) for p in processes]
    for i, key in enumerate(keys):
        print results[i]
        FLTs[key].modelf = np.load(results[i])
        FLTs[key].model = FLTs[key].modelf.reshape((1014,1014))
        os.remove(results[i])
        print FLTs[key].modelf.max()
        
    t1_queue = time.time() 
    #print t1_queue - t0_queue
        
    ######## Slow, test single process execution
    t0_raw = time.time()
    out = _go_compute_model(1, queue=None)
    t1_raw = time.time()
    #print t1_raw - t0_raw
    
    ######## apply_async for getting data out, works the best
    t0_pool = time.time()
    pool = mp.Pool(processes=len(keys))
    results=[pool.apply_async(_go_compute_model, (key, None)) for key in keys]
    pool.close()
    pool.join()
    for res in results:
        key, modelf = res.get(timeout=1)
        FLTs[key].modelf = modelf
        FLTs[key].model = FLTs[key].modelf.reshape((1014,1014))
        print FLTs[key].modelf.max()
        
    t1_pool = time.time()

    print '  Raw: %.3f s (x%d: %.3f)' %(t1_raw - t0_raw, len(keys), (t1_raw - t0_raw)*len(keys))
    print 'Queue: %.3f s (%.1f speedup)' %(t1_queue - t0_queue, (t1_raw - t0_raw)*len(keys)/(t1_queue - t0_queue))
    print ' Pool: %.3f s (%.1f speedup)' %(t1_pool - t0_pool, (t1_raw - t0_raw)*len(keys)/(t1_pool - t0_pool))

### Try pooling with test function looping through object models
def _go_compute_model(ii, queue=None):
    import multiprocessing as mp
    import multiprocessing.queues
    import time

    ### The only way I could figure out getting data in and out and passing through pool.apply_async
    global FLTs
    self = FLTs[ii]
    
    ### Initialize
    self.modelf*=0
    
    ### Loop through objects in given magnitude bins
    #for mag_lim, beams in zip([[10,24], [24,28]], ['ABCDEF', 'A']): 
    for mag_lim, beams in zip([[10,29]], ['A']): 
        ok = (self.catalog['MAG_AUTO'] > mag_lim[0]) & (self.catalog['MAG_AUTO'] < mag_lim[1])
        so = np.argsort(self.catalog['MAG_AUTO'][ok])
        for i in range(ok.sum()):
            ix = so[i]
            #print '%d id=%d mag=%.2f' %(i+1, self.catalog['NUMBER'][ok][ix], self.catalog['MAG_AUTO'][ok][ix])
            for beam in beams:
                self.compute_model(id=self.catalog['NUMBER'][ok][ix], x=self.catalog['x_flt'][ok][ix], y=self.catalog['y_flt'][ok][ix], beam=beam,  sh=[60, 60], verbose=False, in_place=True)
    
    if isinstance(queue, mp.queues.Queue):
        #### For Queue processes
        out = '%s_%.4f.npy' %(mp.current_process().name, time.time())
        np.save(out, self.modelf)
        queue.put(out)
    else:
        #### For pool.apply_async or normal calls
        return (ii, self.modelf)
        
