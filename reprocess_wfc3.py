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

def split_multiaccum(ima, scale_flat=True):
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
            
    NSAMP = ima[0].header['NSAMP']
    sh = ima['SCI',1].shape
    
    cube = np.zeros((NSAMP, sh[0], sh[1]))
    if 'ima' in ima.filename():
        dq = np.zeros((NSAMP, sh[0], sh[1]), dtype=np.int)
    else:
        dq = 0
        
    time = np.zeros(NSAMP)
    for i in range(NSAMP):
        if ima[0].header['UNITCORR'] == 'COMPLETE':
            cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data*ima['TIME',i+1].header['PIXVALUE']/FLAT_F140W
        else:
            cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data/FLAT_F140W
        
        if 'ima' in ima.filename():
            dq[NSAMP-1-i, :, :] = ima['DQ',i+1].data
        
        time[NSAMP-1-i] = ima['TIME',i+1].header['PIXVALUE']
    
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
        
    ### Remove existing products or calwf3 will die
    for ext in ['flt','ima']:
        if os.path.exists(raw.replace('raw', ext)):
            os.remove(raw.replace('raw', ext))
    
    ### Run calwf3
    wfc3tools.calwf3.calwf3(raw)
    
    flt = pyfits.open(raw.replace('raw', 'flt'), mode='update')
    ima = pyfits.open(raw.replace('raw', 'ima'))
    
    cube, dq, time, NSAMP = split_multiaccum(ima, scale_flat=False)
    
    readnoise_2D = np.zeros((1024,1024))
    readnoise_2D[512: ,0:512] += ima[0].header['READNSEA']
    readnoise_2D[0:512,0:512] += ima[0].header['READNSEB']
    readnoise_2D[0:512, 512:] += ima[0].header['READNSEC']
    readnoise_2D[512: , 512:] += ima[0].header['READNSED']
    readnoise_2D = readnoise_2D**2
        
    if len(pop_reads) > 0:
        ### Pop out reads affected by satellite trails or earthshine
        print '\n****\nPop reads %s from %s\n****\n' %(pop_reads, ima.filename())
        
        diff = np.diff(cube, axis=0)
        dt = np.diff(time)
        final_exptime = time[-1]
        final_sci = cube[-1,:,:]*1
        for read in pop_reads:
            final_sci -= diff[read,:,:]
            final_exptime -= dt[read]
        
        #final_var = ima[0].header['READNSEA']**2 + final_sci        
        final_var = readnoise_2D + final_sci        
        final_err = np.sqrt(final_var)/final_exptime
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
        last_ok_read = np.max(zi_flag, axis=0)

        zi_idx = zi < 0
        for i in range(2, NSAMP-1):
            zi_idx[i,:,:] = zi[i,:,:] == last_ok_read

        time_array = time[zi]
        time_array[0,:,:] = 1.e-3 # avoid divide-by-zero
        # pixels that saturated before the last read
        fix = (last_ok_read < (ima[0].header['NSAMP'] - 1)) & (last_ok_read > 1)
        #err = np.sqrt(ima[0].header['READNSEA']**2 + cube)/time_array
        err = np.sqrt(readnoise_2D + cube)/time_array

        final_sci[fix] = np.sum((cube/time_array)*zi_idx, axis=0)[fix]
        final_err[fix] = np.sum(err*zi_idx, axis=0)[fix]

        fixed_sat = (zi_idx.sum(axis=0) > 0) & ((final_dq & 256) > 0)
        final_dq[fixed_sat] -= 256
        print '  Nsat = %d' %(fixed_sat.sum())
        flt['DQ'].data |= final_dq[5:-5,5:-5] & 256
        
    else:
        #### Saturated pixels
        flt['DQ'].data |= ima['DQ',1].data[5:-5,5:-5] & 256
        
    flt['SCI'].data = final_sci[5:-5,5:-5]
    flt['ERR'].data = final_err[5:-5,5:-5]
    
    
    #### Some earthshine flares DQ masked as 32: "unstable pixels"
    mask = (flt['DQ'].data & 32) > 0
    if mask.sum() > 1.e4:
        print '\n****\nTake out excessive DQ=32 flags (N=%e)\n****\n' %(mask.sum())
        flt['DQ'].data[mask] -= 32
        
    ### Update the FLT header
    flt[0].header['IMA2FLT'] = (1, 'FLT extracted from IMA file')
    flt[0].header['IMASAT'] = (fix_saturated*1, 'Manually fixed saturation')

    # if pyfits.__version__ > '3.2':
    #     flt[0].header['IMA2FLT'] = (1, 'FLT extracted from IMA file')
    #     flt[0].header['IMASAT'] = (fix_saturated*1, 'Manually fixed saturation')
    # else:
    #     flt[0].header.update('IMA2FLT', 1, comment='FLT extracted from IMA file')
    #     flt[0].header.update('IMASAT', fix_saturated*1, comment='Manually fixed saturation')
        
    flt.flush()
    
    ### Remove the IMA file
    if remove_ima:
        os.remove(raw.replace('raw', 'ima'))

def show_MultiAccum_reads(raw='ib3701s4q_ima.fits'):
    """
    Make a figure (.ramp.png) showing the individual reads of an 
    IMA or RAW file.
    """    
    img = pyfits.open(raw)
    
    if 'raw' in raw:
        gain=2.5
    else:
        gain=1
        
    cube, dq, time, NSAMP = split_multiaccum(img)
    diff = np.diff(cube, axis=0)
    dt = np.diff(time)
    fig = plt.figure(figsize=[10,10])
    
    #(xs=10, aspect=0.8, wspace=0., hspace=0., left=0.05, NO_GUI=True)
    for j in range(1,NSAMP-1):
        ax = fig.add_subplot(4,4,j)
        ax.imshow(diff[j,:,:]/dt[j]*gain, vmin=0, vmax=4, origin='lower', cmap=plt.get_cmap('hot'))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.text(20,5,'%d' %(j), ha='left', va='bottom', backgroundcolor='white')
    #
    ax = fig.add_subplot(4,4,16)
    ramp_cps = np.median(diff, axis=1)
    avg_ramp = np.median(ramp_cps, axis=1)
    ax.plot(time[2:], ramp_cps[1:,16:-16:4]/100*gain, alpha=0.1, color='black')
    ax.plot(time[2:], avg_ramp[1:]/100*gain, alpha=0.8, color='red', linewidth=2)
    
    fig.tight_layout(h_pad=0.0, w_pad=0.0, pad=0.0)
    plt.savefig(raw.split('_')[0]+'_ramp.png')
    
    return fig
    
    