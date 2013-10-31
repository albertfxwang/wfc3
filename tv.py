"""
Get thermal vac calibration files

Copy grism and imaging flats from /grp/hst/wfc3e/ground_testing/opusdata/ir4
with get_files.sh

cd RAW
grep g1 wfc3_log.prt |grep flat | grep 2008 | sed "s/csii//" |sed "s/_/ /" | sed "s/  /r_/" > list
grep "f[01]" wfc3_log.prt |grep ir | grep flat | grep 2008 | grep spars10 | grep -v spars100| sed "s/csii//" |sed "s/_/ /" |sed "s/  /r_/" >> list


#### Lamp wavelength parameters
dfits *flt.fits |fitsort DATE FILTER EXPTIME OSLAMBDA OSBANDW > lamp_params.info

*** Unfortunately it looks like something was different with teh ground grism configuration compared to 
flight because the order "curtains" don't line up.

"""

import os
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as nd

from threedhst import catIO
import threedhst.dq

from astropy.io import ascii

def testing():
    info = catIO.Readfile("lamp_params.info")
    
    FLAT_F105W = pyfits.open(os.path.join(os.getenv('iref'), 'uc72113oi_pfl.fits'))[1].data[5:-5, 5:-5]
    FLAT_F140W = pyfits.open(os.path.join(os.getenv('iref'), 'uc721143i_pfl.fits'))[1].data[5:-5, 5:-5]
    
    ref = pyfits.open('/Users/brammer/3DHST/Spectra/Work/CONF/sky.G141.set120.fits')
    
    flat = info.filter == 'F140W'
    
    g141 = pyfits.open('RAW/ii140204r_08077120002_flt.fits')
    g141b = pyfits.open('RAW/ii140204r_08074115803_flt.fits')
    
    ds9.view(g141[1].data/FLAT_F140W/1.3e4)
    ds9.view(g141[1].data/FLAT_F105W/1.3e4)
    ds9.view(g141b[1].data/FLAT_F105W/1.3e4)
    
    ds9.view(ref[0].data*FLAT_F140W/FLAT_F105W)
    
def test_regions():
    """
    Make simple model of square "apertures" that are illuminated and shifted for different orders.
    Need +2, +3, +4, +5, 0, -1,-2,-3 orders to get shape right.
    
    dark band at left is from +2 order
    dark band at right is from 0th order
    
    +5th order stronger than others
    
    Note to use "curtain" background in practice should multiply by FLAT_F140W/FLAT_F105W because the bluer 
    flat is more appropriate for the 1.083 um line.
    
    """
    import pyregion
    
    order_scale = [-2.5,-2.5,-2.5,-2.5,-2,-2,-1,-1]
    
    im = pyfits.open('RAW/ii140204r_08077120002_flt.fits')[1]
    orders = np.ones((1014,1014,9))
    s = [5,5,3,3,10,12,6,5] # smooth hard edges of the regions, sigmas in pixels
    scale = np.array([-1,-1,-1,-1,1,1,1,1])*0.03  # slope of spectrum/transmission
    
    #pow = np.ones(9)*0.5
    
    #### Get order masks from region file
    lines = open('orders_sky.reg').readlines()[4:]
    yi, xi = np.indices((1014, 1014))
    for i in range(8):
        fp = open('/tmp/o.reg','w')
        fp.write('image\n%s' %(lines[i]))
        fp.close()
        r = pyregion.open('/tmp/o.reg').as_imagecoord(header=im.header)
        aper = r.get_mask(hdu=im)
        orders[:,:,i] *= nd.filters.gaussian_filter(aper*1., sigma=s[i])*(((xi+1)/507.)*scale[i]+1)
        orders[:,:,i] /= orders[:,:,i].max()
        print i
        
    # M = np.array([[1,1,1,1,0,0,0,0,1],
    #      [1,1,1,1,1,0,0,0,1],
    #      [1,1,1,0,1,0,0,0,1],
    #      [1,1,1,0,1,1,0,0,1],
    #      [1,1,0,0,1,1,0,0,1],
    #      [1,1,0,0,1,1,1,0,1],
    #      [1,0,0,0,1,1,1,0,1],
    #      [1,0,0,0,1,1,1,1,1],
    #      [0,0,0,0,1,1,1,1,1]])
    # 
    # N = np.array([0.89, 1.033, 1.026, 1.038, 1.030, 1.052, 1.050, 1.078, 1.045])
    
    subregion = (ref[0].data*FLAT_F140W/FLAT_F105W)[200:800,:].flatten()
    #subregion = (g141[1].data/FLAT_F105W/1.3e4)[500:800,:].flatten()
    ok = (subregion > 0.5) & (subregion < 1.5)
    M = orders[200:800,:,:].reshape(-1,9)[ok,:]
    N = subregion[ok]
        
    order_scale,res,rank,s = scipy.linalg.lstsq(M, N)
    
    order_model = np.dot(orders, order_scale.reshape(1,9,1)).reshape((1014,1014))
    
    ds9.view(order_model)
    ds9.view(ref[0].data*FLAT_F140W/FLAT_F105W*PAM**-0.2)
    ds9.view(((ref[0].data*FLAT_F140W/FLAT_F105W)-(np.dot(orders, order_scale.reshape(1,9,1)).reshape((1014,1014))*PAM**0.2))+1)
    
    #### Refine +2 (+1 plus 1) order because has different shape
    
    ### Just +2nd order with band at left
    no_p1 = order_scale*1
    no_p1[-2] = 0.
    plus_first = ((ref[0].data*FLAT_F140W/FLAT_F105W)-(np.dot(orders, no_p1.reshape(1,9,1)).reshape((1014,1014))*PAM**0.2))
    
    dy = 15
    
    xarr = np.arange(1014)
    for i in range(0,1000,100):
        plt.plot(xarr+i/900.*dy, np.median(plus_first[i:i+100, :], axis=0), alpha=0.5)
    
    plt.xlim(60,130); plt.ylim(-0.05,0.15)
    
    shift = np.round(xarr/900.*dy-dy/2.)
    sh = plus_first*0.
    for i in range(1014):
        print i
        sh[i,:] += nd.shift(plus_first[i,:], shift[i])
    
    p1_model = np.median(sh, axis=0)
    for i in range(1014):
        orders[i,:,-2] = nd.shift(p1_model, -shift[i])*5
        
    ds9.view(((ref[0].data*FLAT_F140W/FLAT_F105W)-(np.dot(orders, no_p1.reshape(1,9,1)).reshape((1014,1014))*PAM**0.2))+1)
    
def background_components():
    """
    Use high-illum pattern from above for emission line case.
    
    Find a particular observation with pure (high) zodi: COSMOS-11
    """

    curtain = glob.glob('/Users/brammer/3DHST/Spectra/Work/Incoming/ibhj44*flt.fits.gz')
    
    pure_zodi = glob.glob('/Users/brammer/3DHST/Spectra/Work/Incoming/ibhm39*flt.fits.gz')
    
    ims = []
    for im in pure_zodi:
        ims.append(pyfits.open(im))
        
    FLAT_F140W = pyfits.open(os.path.join(os.getenv('iref'), 'uc721143i_pfl.fits'))[1].data[5:-5, 5:-5]
    FLAT_F105W = pyfits.open(os.path.join(os.getenv('iref'), 'uc72113oi_pfl.fits'))[1].data[5:-5, 5:-5]
    FLAT_G141 = pyfits.open(os.path.join(os.getenv('iref'), 'u4m1335mi_pfl.fits'))[1].data[5:-5, 5:-5]
    
    object = ims[0][1].data*FLAT_G141/FLAT_F140W
    dy = 15
    
    xarr = np.arange(1014)
    for i in range(0,1000,100):
        plt.plot(xarr+i/900.*dy, np.median(object[i:i+100, :], axis=0), alpha=0.5)
    #
    shift = np.round(xarr/900.*dy-dy/2.)
    sh = object*0.
    for i in range(1014):
        print i
        sh[i,:] += nd.shift(object[i,:], shift[i])
    #
    model1 = np.median(sh, axis=0)
    model = object*0.
    for i in range(1014):
        model[i,:] = nd.shift(model1, -shift[i])
    
def demo():
    """
    Show dispersed blobs
    """
    curtain = glob.glob('/Users/brammer/3DHST/Spectra/Work/Incoming/ibhj44*flt.fits.gz')
    pure_zodi = glob.glob('/Users/brammer/3DHST/Spectra/Work/Incoming/ibhm39*flt.fits.gz')
    
    ims = []
    for im in curtain[1:2]:
        ims.append(pyfits.open(im))

    for im in pure_zodi[1:2]:
        ims.append(pyfits.open(im))
    
    for i, im in enumerate(ims):
        ds9.frame(i+1)
        ds9.view(im[1].data*FLAT_G141/FLAT_F140W)
        ds9.scale(-1,4)
    