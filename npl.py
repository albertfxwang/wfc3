"""
Noisy power law - Compute statistics of a distribution where the underlying
distribution is a power law that is convolved with Gaussian noise.
"""
import numpy as np

def demo():
    """
    Show test where the source count distribution is a grism model
    """
    import astropy.io.fits as pyfits
    
    flt = pyfits.open('id7h04dvq_GrismFLT.fits')
    flt = pyfits.open('id7h44toq_GrismFLT.fits')
    flt = pyfits.open('id7h01ovq_GrismFLT.fits')
    
    mask = (flt['GDQ'].data == 0) & (flt['GERR'].data < 0.1) & (flt['GERR'].data > 0)
    
    noise = np.random.normal(size=mask.sum())*np.median(flt['GERR'].data[mask])
    data = flt['MODEL'].data[mask] + noise
    
    R = np.arange(-10,10,1.e-4)
    noise = np.random.normal(size=R.size)*ss
    data = np.exp(-R**2/2/sig**2)**2 + noise
    
    from astropy.modeling import models
    n = 1
    r_e = 100
    ser = models.Sersic2D(amplitude=1, r_eff=r_e, n=n, x_0=500, y_0=500, ellip=0, theta=0)
    yp, xp = np.indices((1014,1014))
    ser_data = ser(xp, yp)
    
    ser_data /= ser_data.max()
    
    noise = np.random.normal(size=ser_data.shape)*ss
    
    data = (ser_data + noise).flatten()
    
    #data *= 100
    #data = noise
    
    import mywfc3.npl 
    self = mywfc3.npl.NoisyPowerlaw(data=data, c0=[0,0,0,0,-2])
    self.dx_min = 0
    x = np.arange(self.range[0], self.range[1], self.s0/50)
    res, fit = self.fit(x=x, tol=1.e-6)
    plt.plot(self.xh, self.yh, color='k', linestyle='steps-mid')
    plt.plot(self.xh, fit, color='r', linestyle='steps-mid')
    
    sn_mask =  (flt['MODEL'].data < 0.5*flt['GERR'].data)[mask]
    
class NoisyPowerlaw(object):
    def __init__(self, data=[], range=None, m0=None, s0=None, c0=[0,-2], log=False, hist_bins='knuth', sigma_range=[-4,40], flt_file=None, flt_ext=1):
        """"Statistics of a noisy powerlaw distribution
        
        Compute statistics of a distribution where the underlying distribution 
        is a power law that is convolved with Gaussian noise.
        
        Parameters
        ----------
        data : type
        
        range : type
        
        
        Returns
        -------

        """
        import astropy.stats
        import astropy.io.fits as pyfits
        
        ## Get data from an FLT file
        if flt_file is not None:
            im = pyfits.open(flt_file)

            if '_flc' in im.filename():
                imscl = im[0].header['EXPTIME']
            else:
                imscl = 1
            
            # DQ mask and reasonable uncertainties
            mask = ((im['DQ',flt_ext].data == 0) & 
                    (im['ERR',flt_ext].data < 1*imscl))
                    
            data = im['SCI',flt_ext].data[mask]
            if mask.sum() > 1.e6:
                rnd_ix = np.argsort(np.random.random(size=mask.sum()))
                data = data[rnd_ix[:int(1e6)]]
            
        self.data = data
        
        if m0 is None:
            self.m0 = np.median(data)
        else:
            self.m0 = m0
            
        if s0 is None:
            self.s0 = astropy.stats.mad_std(data)
        else:
            self.s0 = s0
        
        self.dx_min = self.s0/10
        
        self.c0 = c0
        
        if range is None:
            self.range = [self.m0+sigma_range[0]*self.s0, 
                          self.m0+sigma_range[1]*self.s0]
        else:
            self.range = range
        
        self.log = log
        if log:
            h = astropy.stats.histogram(np.log(data), bins=hist_bins, range=np.log(self.range), normed=False)
            self.bins = h[1]
            self.xh = np.exp(h[1][1:] - np.diff(h[1])/2.)
            self.yh = h[0]
        else:
            h = astropy.stats.histogram(data, bins=hist_bins, range=self.range, normed=False)
            self.bins = h[1]
            self.xh = h[1][1:] - np.diff(h[1])/2.
            self.yh = h[0]
        
        self.sum = self.yh.sum()
        self.bin_width = np.diff(self.bins)
        
        if False:
            self = npl.NoisyPowerlaw(data, log=False)
            
            x = np.arange(self.range[0], self.range[1], self.s0/10)
            #x = self.xh
            dx = np.diff(x)[0]
            
            params = [-2*10, 0, 0, self.m0-0.015, self.s0/dx*0.95]
            smx = self._objective(params, 10, x, self, 1, True)
            res, model = self.fit(x=x, guess=params, scl=10)
            
            plt.plot(self.xh, self.yh, color='k', linewidth=1, alpha=0.5, linestyle='steps-mid')
            plt.plot(self.xh, model, color='r', linewidth=1, alpha=0.5, linestyle='steps-mid')
            
            #smf = self._objective(res, 10, x, self, 1)
            
    def fit(self, x=None, guess=None, method='powell', scl=10, verbose=True, **kwargs):
        """TBD
        
        Parameters
        ----------
        guess : type
        
        scl : type
            
        verbose : bool
        
        **kwargs passed to `~scipy.optimize.fmin_powell`.
        
        Returns
        -------


        """
        if x is None:
            x = np.arange(self.range[0], self.range[1], self.s0/50)
                        
        if guess is None:
            guess = np.hstack((np.array(self.c0), self.m0-0.015, self.s0*0.95))
        
        Np = len(guess)-2
        power = np.hstack((scl**np.arange(1,Np+1)[::-1], scl**2, scl**2))
        init = np.array(guess)*power
        
        import scipy.optimize
        res = scipy.optimize.minimize(self._objective, init, args=(scl, x, self, 0, verbose), method=method, **kwargs)
        
        model = self._objective(res.x, scl, x, self, 1, True)
        res.x /= power
        
        #rms = res.x[-1]*np.diff(x[:2])[0]
        return res, model
        
           
    @staticmethod
    def _objective(params, scl, x, self, ret, verbose):
        """Objective function for the fit
        
        Parameters
        ----------
        params : type
        
        scl : type
        
        x : type
        
        self : class
        
        ret : int
        
        verbose : bool
        Returns
        -------
        If `ret==0`, return chi-squared, else return the fit function.
        
        """
        
        import scipy.ndimage as nd
        
        #power = scl**np.arange(1,len(params)-1)
        Np = len(params)-2
        power = np.hstack((scl**np.arange(1,Np+1)[::-1], scl**2, scl**2))
        
        coeffs = np.array(params[:-2])/power[:-2]
        x0 = params[-2]/power[-2]
        y = NoisyPowerlaw.shifted_powerlaw(x, x0=x0, coeffs=coeffs, dx_min=self.dx_min)
        
        dx = np.diff(x[:2])[0]
        
        s = params[-1]/power[-1]
        #print(s/dx)
        if s < 0:
            s = 0.01
            
        sm = nd.gaussian_filter1d(y, s/dx)
        
        # Interpolate to histogram bins
        sm = np.interp(self.xh, x, sm)
        sm *= self.sum*self.bin_width/np.trapz(sm, self.xh)
        
        err = np.sqrt(self.yh + (0.02*self.yh)**2)
        
        ok = self.yh > 0
        chi2 = np.sum((sm-self.yh)[ok]**2/err[ok]**2)
        
        if verbose:
            print(params/power, chi2)
        
        if ret == 0:
            return chi2
        else:
            return sm
        
    @classmethod
    def shifted_powerlaw(cls, x, x0=0.9, coeffs=[0, -2], dx_min=0.01):
        """Shifted powerlaw function
        
        Return a shifted powerlaw function with exponential coefficients 
        `coeffs`:
        
            >>> lnx = np.log(x-x0) # where x > x0
            >>> lny = polyval(coeffs, lnx)
            >>> y = np.exp(lny) # and set to zero where x <= x0
        
        Parameters
        ----------
        x : `~numpy.ndarray`
            Array on which to compute the function.
            
        x0 : float
            "zero" of the distribution
            
        coeffs : list of floats
            Powerlaw coefficients, in reverse order as for `~scipy.polyval`:
        
        dx_min : float
            Limit to control exponential behavior of the powerlaw as x -> x0.
            Clip values at (x0, x0+dx_min) to y[x0+dx_min]:
            
                >>> y[(x > x0) & (x < (x0+dx_min))] = y[x >= (x0+dx_min)][0]
                
        Returns
        -------
        y : `~numpy.ndarray`  
        """
        from scipy import polyval
        
        y = x*0.

        ok = x > x0
        lnx = np.log(x[ok]-x0)
        lny = polyval(np.hstack((coeffs, 0)), lnx)
        y[ok] = np.exp(lny)
        
        clip = (x > x0) & (x < (x0+dx_min))
        y[clip] = y[x >= (x0+dx_min)][0] #/clip.sum()

        return y
     
            
            