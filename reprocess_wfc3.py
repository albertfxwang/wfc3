"""
Scripts to reprocess WFC3 exposures with time-variable backgrounds 
or satellite trails.  
"""

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np
import matplotlib.pyplot as plt

import os
import glob
import shutil

import logging
logger = logging.getLogger('reprocess_wfc3')
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(message)s', datefmt='%m/%d/%Y %H:%M:%S')
ch.setFormatter(formatter)
if len(logger.handlers) == 0:
    logger.addHandler(ch)

#logger.info('test')

#logging.basicConfig(format='%(name)s, %(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

def split_multiaccum(ima, scale_flat=True, get_err=False):
    """
    Pull out the MultiAccum reads of a RAW or IMA file into a single 3D 
    matrix.
    
    Returns cube[NSAMP,1024,1014], time, NSAMP
    """    
    skip_ima = ('ima' in ima.filename()) & (ima[0].header['FLATCORR'] == 'COMPLETE')
    if scale_flat & ~skip_ima:
        FLAT_F140W = pyfits.open(os.path.join(os.getenv('iref'), 'uc721143i_pfl.fits'))[1].data
    else:
        FLAT_F140W = 1
    
    
    is_dark = 'drk' in ima.filename()
    if is_dark:
        FLAT_F140W = 1
                
    NSAMP = ima[0].header['NSAMP']
    sh = ima['SCI',1].shape
    
    cube = np.zeros((NSAMP, sh[0], sh[1]))
    if 'ima' in ima.filename():
        dq = np.zeros((NSAMP, sh[0], sh[1]), dtype=np.int)
    else:
        dq = 0
    
    if get_err:
        cube_err = cube*0
        
    time = np.zeros(NSAMP)
    for i in range(NSAMP):
        if (ima[0].header['UNITCORR'] == 'COMPLETE') & (~is_dark):
            cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data*ima['TIME',i+1].header['PIXVALUE']/FLAT_F140W
        else:
            #print 'Dark'
            cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data/FLAT_F140W
        
        if get_err:
            if ima[0].header['UNITCORR'] == 'COMPLETE':
                cube_err[NSAMP-1-i, :, :] = ima['ERR',i+1].data*ima['TIME',i+1].header['PIXVALUE']/FLAT_F140W
            else:
                cube_err[NSAMP-1-i, :, :] = ima['ERR',i+1].data/FLAT_F140W
            
        if 'ima' in ima.filename():
            dq[NSAMP-1-i, :, :] = ima['DQ',i+1].data
        
        time[NSAMP-1-i] = ima['TIME',i+1].header['PIXVALUE']
    
    if get_err:
        return cube, cube_err, dq, time, NSAMP
    else:
        return cube, dq, time, NSAMP
        
def make_IMA_FLT(raw='ibhj31grq_raw.fits', pop_reads=[], remove_ima=True, fix_saturated=True):
    """
    Run calwf3, if necessary, to generate ima & flt files.  Then put the last
    read of the ima in the FLT SCI extension and let Multidrizzle flag 
    the CRs.
    
    Optionally pop out reads affected by satellite trails or earthshine.  The 
    parameter `pop_reads` is a `list` containing the reads to remove, where
    a value of 1 corresponds to the first real read after the 2.9s flush.
    
    Requires IRAFX for wfc3tools
    """
    import wfc3tools
        
    #### Remove existing products or calwf3 will die
    for ext in ['flt','ima']:
        if os.path.exists(raw.replace('raw', ext)):
            os.remove(raw.replace('raw', ext))
    
    #### Run calwf3
    wfc3tools.calwf3.calwf3(raw)
    
    flt = pyfits.open(raw.replace('raw', 'flt'), mode='update')
    ima = pyfits.open(raw.replace('raw', 'ima'))
    
    #### Pull out the data cube, order in the more natural sense
    #### of first reads first
    cube, dq, time, NSAMP = split_multiaccum(ima, scale_flat=False)
    
    #### Readnoise in 4 amps
    readnoise_2D = np.zeros((1024,1024))
    readnoise_2D[512: ,0:512] += ima[0].header['READNSEA']
    readnoise_2D[0:512,0:512] += ima[0].header['READNSEB']
    readnoise_2D[0:512, 512:] += ima[0].header['READNSEC']
    readnoise_2D[512: , 512:] += ima[0].header['READNSED']
    readnoise_2D = readnoise_2D**2

    #### Gain in 4 amps
    gain_2D = np.zeros((1024,1024))
    gain_2D[512: ,0:512] += ima[0].header['ATODGNA']
    gain_2D[0:512,0:512] += ima[0].header['ATODGNB']
    gain_2D[0:512, 512:] += ima[0].header['ATODGNC']
    gain_2D[512: , 512:] += ima[0].header['ATODGND']
    
    ### Pop out reads affected by satellite trails or earthshine
    if len(pop_reads) > 0:
        print '\n****\nPop reads %s from %s\n****\n' %(pop_reads, ima.filename())
        
        #### Need to put dark back in for Poisson
        dark_file = ima[0].header['DARKFILE'].replace('iref$', os.getenv('iref')+'/')
        dark = pyfits.open(dark_file)
        dark_cube, dark_dq, dark_time, dark_NSAMP = split_multiaccum(dark, scale_flat=False)
        
        #### Need flat for Poisson
        flat_file = ima[0].header['PFLTFILE'].replace('iref$', os.getenv('iref')+'/')
        flat = pyfits.open(flat_file)#[1].data
        ff = flat[1].data
        
        #### Subtract diffs if flagged reads
        diff = np.diff(cube, axis=0)
        dark_diff = np.diff(dark_cube, axis=0)

        dt = np.diff(time)
        final_exptime = time[-1]
        final_sci = cube[-1,:,:]*1
        final_dark = dark_cube[NSAMP-1,:,:]*1        
        for read in pop_reads:
            final_sci -= diff[read,:,:]
            final_dark -= dark_diff[read,:,:]
            final_exptime -= dt[read]
                
        #### Variance terms
        ## read noise
        final_var = readnoise_2D*1
        ## poisson term
        final_var += (final_sci*ff + final_dark*gain_2D)*(gain_2D/2.368)
        ## flat errors
        final_var += (final_sci*ff*flat['ERR'].data)**2
        final_err = np.sqrt(final_var)/ff/(gain_2D/2.368)/1.003448/final_exptime
        
        final_sci /= final_exptime
                
        flt[0].header['EXPTIME'] = final_exptime
        
    else:
        final_sci = ima['SCI', 1].data*1
        final_err = ima['ERR', 1].data*1
    
    final_dq = ima['DQ', 1].data*1
    
    #### For saturated pixels, look for last read that was unsaturated
    #### Background will be different under saturated pixels but maybe won't
    #### matter so much for such bright objects.
    if fix_saturated:
        print 'Fix Saturated pixels:'
        #### Saturated pixels
        zi, yi, xi = np.indices(dq.shape)
        saturated = (dq & 256) > 0
        # 1024x1024 index array of reads where pixels not saturated
        zi_flag = zi*1
        zi_flag[saturated] = 0
        ### 2D array of the last un-saturated read
        last_ok_read = np.max(zi_flag, axis=0)
        sat_zero = last_ok_read == 0        
        pyfits.writeto(raw.replace('_raw','_lastread'), data=last_ok_read[5:-5,5:-5], header=flt[1].header, clobber=True)
        ### keep pixels from first read even if saturated
        last_ok_read[sat_zero] = 1
        
        zi_idx = zi < 0
        for i in range(1, NSAMP-1):
            zi_idx[i,:,:] = zi[i,:,:] == last_ok_read

        time_array = time[zi]
        time_array[0,:,:] = 1.e-3 # avoid divide-by-zero
        # pixels that saturated before the last read
        fix = (last_ok_read < (ima[0].header['NSAMP'] - 1)) & (last_ok_read > 0)
        #err = np.sqrt(ima[0].header['READNSEA']**2 + cube)/time_array
        err = np.sqrt(readnoise_2D + cube)/time_array

        final_sci[fix] = np.sum((cube/time_array)*zi_idx, axis=0)[fix]
        final_err[fix] = np.sum(err*zi_idx, axis=0)[fix]

        fixed_sat = (zi_idx.sum(axis=0) > 0) & ((final_dq & 256) > 0)
        final_dq[fixed_sat] -= 256
        final_dq[sat_zero] |= 256
        
        print '  Nsat = %d' %(fixed_sat.sum())
        flt['DQ'].data |= final_dq[5:-5,5:-5] & 256
        
    else:
        #### Saturated pixels
        flt['DQ'].data |= ima['DQ',1].data[5:-5,5:-5] & 256
        
    flt['SCI'].data = final_sci[5:-5,5:-5]
    flt['ERR'].data = final_err[5:-5,5:-5]
    
    #### Some earthshine flares DQ masked as 32: "unstable pixels"
    mask = (flt['DQ'].data & 32) > 0
    if mask.sum() > 2.e4:
        print '\n****\nTake out excessive DQ=32 flags (N=%e)\n****\n' %(mask.sum())
        flt['DQ'].data[mask] -= 32
        
    ### Update the FLT header
    flt[0].header['IMA2FLT'] = (1, 'FLT extracted from IMA file')
    flt[0].header['IMASAT'] = (fix_saturated*1, 'Manually fixed saturation')
    flt[0].header['NPOP'] = (len(pop_reads), 'Number of reads popped from the sequence')
    for iread, read in enumerate(pop_reads):
        flt[0].header['POPREAD%d' %(iread+1)] = (read, 'Read kicked out of the MULTIACCUM sequence')
        
    flt.flush()
    
    ### Remove the IMA file
    if remove_ima:
        os.remove(raw.replace('raw', 'ima'))

def show_MultiAccum_reads(raw='ibp329isq_raw.fits', flatten_ramp=False, verbose=True):
    """
    Make a figure (.ramp.png) showing the individual reads of an 
    IMA or RAW file.
    """    
    import scipy.ndimage as nd
    
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)
        
    img = pyfits.open(raw)
    
    if 'raw' in raw:
        gains = [2.3399999, 2.3699999, 2.3099999, 2.3800001]
        gain = np.zeros((1024,1024))
        gain[512: ,0:512] += gains[0]
        gain[0:512,0:512] += gains[1]
        gain[0:512, 512:] += gains[2]
        gain[512: , 512:] += gains[3]
    else:
        gain=1
    
    logger.info('Make MULTIACCUM cube')
        
    #### Split the multiaccum file into individual reads    
    cube, dq, time, NSAMP = split_multiaccum(img, scale_flat=False)
    
    if 'raw' in raw:
        dark_file = img[0].header['DARKFILE'].replace('iref$', os.getenv('iref')+'/')
        dark = pyfits.open(dark_file)
        dark_cube, dark_dq, dark_time, dark_NSAMP = split_multiaccum(dark, scale_flat=False)

        diff = np.diff(cube-dark_cube[:NSAMP,:,:], axis=0)*gain
        dt = np.diff(time)
    
        #### Need flat for Poisson
        flat_file = img[0].header['PFLTFILE'].replace('iref$', os.getenv('iref')+'/')
        flat = pyfits.open(flat_file)#[1].data
        ff = flat[1].data
        diff /= ff
    else:
        diff = np.diff(cube, axis=0)
        dt = np.diff(time)
    
    ####  Average ramp
    ramp_cps = np.median(diff, axis=1)
    avg_ramp = np.median(ramp_cps, axis=1)
    
    #### Initialize the figure
    logger.info('Make plot')
    
    plt.ioff()
    fig = plt.figure(figsize=[10,10])

    ## Smoothing
    smooth = 1
    kernel = np.ones((smooth,smooth))/smooth**2
    
    ## Plot the individual reads
    for j in range(1,NSAMP-1):
        ax = fig.add_subplot(4,4,j)
        smooth_read = nd.convolve(diff[j,:,:],kernel)
        ax.imshow(smooth_read[5:-5:smooth, 5:-5:smooth]/dt[j], 
                  vmin=0, vmax=4, origin='lower', cmap=plt.get_cmap('hot'))
        
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(20,5,'%d' %(j), ha='left', va='bottom', backgroundcolor='white')
    
    ## Show the ramp
    ax = fig.add_axes((0.6, 0.05, 0.37, 0.18))
    ax.plot(time[2:], (ramp_cps[1:,16:-16:4].T/np.diff(time)[1:]).T, alpha=0.1, color='black')
    ax.plot(time[2:], avg_ramp[1:]/np.diff(time)[1:], alpha=0.8, color='red', linewidth=2)
    ax.set_xlabel('time'); ax.set_ylabel('background [e/s]')

    fig.tight_layout(h_pad=0.3, w_pad=0.3, pad=0.5)
    root=raw.split('_')[0]
    plt.savefig(root+'_ramp.png')
    
    #### Same ramp data file    
    np.savetxt('%s_ramp.dat' %(root), np.array([time[1:], avg_ramp/np.diff(time)]).T, fmt='%.3f')
    
    if flatten_ramp:
        #### Flatten the ramp by setting background countrate to the average.  
        #### Output saved to "*x_flt.fits" rather than the usual *q_flt.fits.
        import wfc3tools
        
        flux = avg_ramp/np.diff(time)
        avg = avg_ramp.sum()/time[-1]
        min = flux[1:].min()
        subval = np.cumsum((flux-avg)*np.diff(time))
        
        imraw = pyfits.open(raw.replace('ima','raw'))
        for i in range(1, NSAMP):
            logger.info('Remove excess %.2f e/s from read #%d (t=%.1f)' %(flux[-i]-min, NSAMP-i+1, time[-i]))
            
            imraw['SCI',i].data = imraw['SCI',i].data - np.cast[int](subval[-i]/2.36*ff)
                
        files=glob.glob(raw.split('q_')[0]+'x_*')
        for file in files:
            os.remove(file)
            
        imraw[0].header['CRCORR'] = 'PERFORM'
        imraw.writeto(raw.replace('q_raw', 'x_raw'), clobber=True)
        
        ## Run calwf3
        wfc3tools.calwf3.calwf3(raw.replace('q_raw', 'x_raw'))
                
    return fig
    
def clean_ramp(raw='ibp329isq_raw.fits'):
    tr, fr = np.loadtxt(raw.replace('raw.fits','ramp.dat'), unpack=True)
    
def ramp_crrej(raw):
    """
    testing xxx
    
    Do own CR-rejection
    """
    import wfc3tools
    
    ### Run calwf3
    if not os.path.exists(raw.replace('raw', 'ima')):
        ### Remove existing products or calwf3 will die
        if os.path.exists(raw.replace('raw', 'flt')):
            os.remove(raw.replace('raw', 'flt'))
        
        wfc3tools.calwf3.calwf3(raw)
        
    ima = pyfits.open(raw.replace('raw', 'ima'))
    
    ### cube: cumulated reads
    cube, cube_err, dq, time, NSAMP = split_multiaccum(ima, scale_flat=False, get_err=True)
    ### diff: Delta reads
    diff = np.diff(cube, axis=0)
    dt = np.diff(time)
    
    ### units right? should be e-/s
    err = np.sqrt((cube_err[1:,:,:].T/time[1:])**2 + (cube_err[:-1,:,:].T/time[:-1])**2).T
    # ix = -8
    # med = np.median((diff[ix,400:600,400:600]/dt[ix]))
    # rnd = np.random.normal(size=(1024,1024))*err[ix,:,:]+med
    
    #### Compute average ramp
    ## bit mask
    mask = (dq[-1,:,:] & ~(576+8192+4096+32)) == 0
    ## bright objects
    mask &= np.abs((diff[-1,:,:]-np.median(diff[-1,:,:]))/dt[-1]) < 5*err[-1,:,:]
    
    #### countrate:  delta countrate, minus average ramp
    countrate = diff*0.
    ramp = dt*0.
    for i in range(len(dt)):
        ramp[i] = np.median(diff[i,:,:][mask])/dt[i]
        countrate[i,:,:] = diff[i,:,:]/dt[i]-ramp[i]
    
    #### CR rejection:  (countrate_i - med)/err > threshold    
    med = np.median(countrate, axis=0)
    cr_full = ~np.isfinite(med)
    sum = cr*0.
    sum_time = sum*0.
    for i in range(5,len(dt)):
        cr = (countrate[i,:,:]-med)/err[i,:,:] < 15
        cr_full |= cr
        sum += countrate[i,:,:]*dt[i]*cr
        sum_time += dt[i]*cr
    
def split_ima():
    
    import scipy.ndimage as nd
    
    ima_file = 'ibp329isq_ima.fits'
    ima = pyfits.open(ima_file)   
     
    #### Split the multiaccum file into individual reads    
    cube, cube_err, dq, time, NSAMP = mywfc3.reprocess_wfc3.split_multiaccum(ima, scale_flat=False, get_err=True)
        
    #
    dark_file = ima[0].header['DARKFILE'].replace('iref$', os.getenv('iref')+'/')
    dark = pyfits.open(dark_file)
    dark_cube, dark_dq, dark_time, dark_NSAMP = mywfc3.reprocess_wfc3.split_multiaccum(dark, scale_flat=False)
    
    #
    #### Readnoise in 4 amps
    readnoise_2D = np.zeros((1024,1024))
    readnoise_2D[512: ,0:512] += ima[0].header['READNSEA']
    readnoise_2D[0:512,0:512] += ima[0].header['READNSEB']
    readnoise_2D[0:512, 512:] += ima[0].header['READNSEC']
    readnoise_2D[512: , 512:] += ima[0].header['READNSED']
    readnoise_2D = readnoise_2D**2

    #### Gain in 4 amps
    gain_2D = np.zeros((1024,1024))
    gain_2D[512: ,0:512] += ima[0].header['ATODGNA']
    gain_2D[0:512,0:512] += ima[0].header['ATODGNB']
    gain_2D[0:512, 512:] += ima[0].header['ATODGNC']
    gain_2D[512: , 512:] += ima[0].header['ATODGND']
    
    #### Need flat for Poisson
    flat_file = ima[0].header['PFLTFILE'].replace('iref$', os.getenv('iref')+'/')
    flat = pyfits.open(flat_file)#[1].data
    ff = flat[1].data
    
    #### Subtract diffs if flagged reads
    diff = np.diff(cube, axis=0)
    dark_diff = np.diff(dark_cube, axis=0)
    
    dt = np.diff(time)
    
    diff_var = diff*0.
    diff_err = diff*0.
    for i in range(diff.shape[0]):
        diff_var[i,:,:] = readnoise_2D*1.
        diff_var[i,:,:] += (diff[i,:,:]*ff + dark_diff[i,:,:]*gain_2D)*(gain_2D/2.368)
        diff_var[i,:,:] += (diff[i,:,:]*ff*flat['ERR'].data)**2
        diff_err[i,:,:] = np.sqrt(diff_var[i,:,:])/ff/(gain_2D/2.368)/1.003448
    
    flt = pyfits.open(ima_file.replace('ima','flt'))
    
    letters = 'abcdefghijklmnoprstuv'
    for i in range(1, diff.shape[0]):
        dq_i = dq[i+1]
        crset = (dq_i & 8192) > 0
        dq_i[crset] -= 8192
        flt['SCI'].data = diff[i,5:-5,5:-5]/dt[i]
        flt['ERR'].data = diff_err[i,5:-5,5:-5]/dt[i]
        flt['DQ'].data = dq_i[5:-5,5:-5]
        flt.writeto(flt.filename().replace('q_', letters[i-1]+'_'), clobber=True)
        
    import drizzlepac
    read_files=glob.glob('*[a-p]_flt.fits')
    drizzlepac.astrodrizzle.AstroDrizzle(read_files, output='test', clean=True, context=False, preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, driz_cr_corr=True, driz_combine=True, resetbits=0)
    
    drizzlepac.astrodrizzle.AstroDrizzle(flt.filename(), output='testx', clean=True, context=False, preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, driz_cr_corr=True, driz_combine=True, resetbits=0)
    
    
        