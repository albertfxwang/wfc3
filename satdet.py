"""
Satellite trail detection, inspired by Borncamp & Lam ACS-ISR 16-01

For WFC3/IR, use rather the DQ/time arrays to identified satellite 
trails flagged as CRs but with masks that need to be expanded.
"""

# STDLIB
import glob
import multiprocessing
import time
import warnings
from multiprocessing import Process, Queue

# THIRD PARTY
import numpy as np
from astropy.io import fits
#from astropy.stats import biweight_location
from astropy.stats import biweight_midvariance, sigma_clipped_stats
from astropy.utils.exceptions import AstropyUserWarning

import skimage
from astropy.utils.introspection import minversion
from scipy import stats, polyfit, polyval
from skimage import transform
from skimage import morphology as morph
from skimage.measure import label, regionprops

from skimage import exposure
from skimage.feature import canny

def _detsat_one(filename, ext=1, small_edge=12, buf=8, min_len=100, min_e=0.8, min_area=500, plot=False, verbose=False, ds9=None, update=0):
    """Called by :func:`detsat`."""
    if verbose:
        t_beg = time.time()

    fname = '{0}[{1}]'.format(filename, ext)

    # check extension
    if ext not in [1, ('SCI',1)]: #in (1, 4, 'SCI', ('SCI', 1), ('SCI', 2)):
        warnings.warn('{0} is not a valid science extension for '
                      'ACS/WFC'.format(ext), AstropyUserWarning)

    # get the data
    im = fits.open(filename)
    image = im[ext].data
    
    imtime = im['TIME'].data
    dq = im['DQ'].data
    
    edge = (imtime < imtime.max()) & (dq == 0)        
    morph.remove_small_objects(edge, min_size=small_edge, 
                               connectivity=small_edge*3,
                               in_place=True)  
                               
    edge = morph.binary_dilation(edge, selem = np.ones((buf,buf)))
    #edge = morph.binary_erosion(edge, selem = np.ones((3,3)))
    
    if ds9 is not None:
        ds9.frame(1); ds9.view(edge*1)
    
    # Label the image and remove short, round features    
    label_image = label(edge)
    regs = regionprops(label_image)
    if verbose:
        print('label area e length')
    
    for reg in regs:
        #if verbose:
        #    print('{0:3d} {1:6d} {2:0.2f} {3:0.2f}'.format(reg.label, reg.area, reg.eccentricity, reg.major_axis_length))
            
        if (reg.area < min_area) | (reg.eccentricity < min_e) | (reg.major_axis_length < min_len):
            for coo in reg.coords:
                label_image[coo[0], coo[1]] -= reg.label
    
    edge = label_image > 0
    
    ## Fit lines to cleaned labels
    label_image = label(edge)
    regs = regionprops(label_image)
    yp, xp = np.indices(imtime.shape)
    ntrail = 0
    for reg in regs:
        if reg.eccentricity < 0.95:
            continue
            
        c = polyfit(reg.coords[:,1], reg.coords[:,0], 1)
        
        yline = polyval(c, xp)
        theta = np.arctan(c[0])
        ywidth = reg.minor_axis_length/np.cos(theta)*0.8
        line_mask = np.abs(yp-yline) < ywidth
        if verbose:
            print('{0:3d} {1:6d} {4:6d} {2:0.2f} {3:0.2f}'.format(reg.label, reg.area, reg.eccentricity, reg.major_axis_length, line_mask.sum()))
        
        if line_mask.sum() > 10*reg.area:
            continue
        else:
            ntrail += 1
                
        edge |= line_mask
    
    if ds9 is not None:
        ds9.frame(2); ds9.view(edge*1)
        ds9.frame(3); ds9.view((image-np.median(image[dq == 0]))*((~edge) & (dq == 0)))
    
    if update:
        im = fits.open(filename, mode='update')
        im[0].header['SATTRAIL'] = (ntrail, 'Number of detected satellite trails')
        
        im['DQ'].data |= update*edge
        
        im.flush()
        
    return edge
    