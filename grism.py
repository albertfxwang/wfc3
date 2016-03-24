"""
Demonstrate aXe trace polynomials.

  v1.0 - October 14, 2014  (G. Brammer, N. Pirzkal, R. Ryan) 

"""
import numpy as np

class aXeConf():
    def __init__(self, conf_file='WFC3.IR.G141.V2.5.conf'):
        
        if conf_file is not None:
            self.conf = self.read_conf_file(conf_file)
            self.conf_file = conf_file
            self.count_beam_orders()
            
            ## Global XOFF/YOFF offsets
            if 'XOFF' in self.conf.keys():
                self.xoff = np.float(conf['XOFF'])
            else:
                self.xoff = 0.

            if 'YOFF' in self.conf.keys():
                self.yoff = np.float(conf['YOFF'])
            else:
                self.yoff = 0.
            
    def read_conf_file(self, conf_file='WFC3.IR.G141.V2.5.conf'):
        """
        Read an aXe config file, convert floats and arrays
        """
        from collections import OrderedDict
    
        conf = OrderedDict()
        lines = open(conf_file).readlines()
        for line in lines:
            ## empty / commented lines
            if (line.startswith('#')) | (line.strip() == '') | ('"' in line):
                continue
        
            ## split the line, taking out ; and # comments
            spl = line.split(';')[0].split('#')[0].split()
            param = spl[0]
            if len(spl) > 2: 
                value = np.cast[float](spl[1:])
            else:
                try:
                    value = float(spl[1])
                except:
                    value = spl[1]

            conf[param] = value
    
        return conf
    
    def count_beam_orders(self):
        """
        Get the maximum polynomial order in DYDX or DLDP for each beam
        """
        self.orders = {}
        for beam in ['A','B','C','D','E','F','G','H','I','J']:
            order = 0
            while 'DYDX_%s_%d' %(beam, order) in self.conf.keys():
                order += 1
            
            while 'DLDP_%s_%d' %(beam, order) in self.conf.keys():
                order += 1
            
            self.orders[beam] = order-1

    def get_beams(self):
        """
        Get beam parameters and sensitivity curves
        """
        import os
        import collections
        from astropy.table import Table
        
        self.dxlam = collections.OrderedDict()
        self.nx = collections.OrderedDict()
        self.sens = collections.OrderedDict()
        self.beams = []
        
        for beam in self.orders:
            if self.orders[beam] > 0:
                self.beams.append(beam)
                self.dxlam[beam] = np.arange(self.conf['BEAM%s' %(beam)][0], self.conf['BEAM%s' %(beam)][1], dtype=int)
                self.nx[beam] = int(self.dxlam[beam].max()-self.dxlam[beam].min())+1
                self.sens[beam] = Table.read('%s/%s' %(os.path.dirname(self.conf_file), self.conf['SENSITIVITY_%s' %(beam)]))
                self.sens[beam].wave = np.cast[np.double](self.sens[beam]['WAVELENGTH'])
                self.sens[beam].sens = np.cast[np.double](self.sens[beam]['SENSITIVITY'])
    
                
    def field_dependent(self, xi, yi, coeffs):
        """
        aXe field-dependent coefficients
        """
        ## number of coefficients for a given polynomial order
        ## 1:1, 2:3, 3:6, 4:10, order:order*(order+1)/2
        if isinstance(coeffs, float):
            order = 1
        else:
            order = int(-1+np.sqrt(1+8*len(coeffs)))/2
    
        ## Build polynomial terms array
        ## $a = a_0+a_1x_i+a_2y_i+a_3x_i^2+a_4x_iy_i+a_5yi^2+$ ...
        xy = []
        for p in range(order):
            for px in range(p+1):
                xy.append(xi**(p-px)*yi**(px))
    
        ## Evaluate the polynomial, allowing for N-dimensional inputs
        a = np.sum((np.array(xy).T*coeffs).T, axis=0)
    
        return a
    
    def get_beam_trace(self, x=507, y=507, dx=0., beam='A'):
        """
        Get an aXe beam trace for an input reference pixel and 
        list of output x pixels dx
        """
        NORDER = self.orders[beam]+1
        
        xi, yi = x-self.xoff, y-self.yoff
        xoff_beam = self.field_dependent(xi, yi, self.conf['XOFF_%s' %(beam)])
        yoff_beam = self.field_dependent(xi, yi, self.conf['YOFF_%s' %(beam)])
    
        ## y offset of trace (DYDX)
        dydx = np.zeros(NORDER) #0 #+1.e-80
        for i in range(NORDER):
            if 'DYDX_%s_%d' %(beam, i) in self.conf.keys():
                coeffs = self.conf['DYDX_%s_%d' %(beam, i)]
                dydx[i] = self.field_dependent(xi, yi, coeffs)
            
        # $dy = dydx_0+dydx_1 dx+dydx_2 dx^2+$ ...
        dy = yoff_beam
        for i in range(NORDER):
            dy += dydx[i]*(dx-xoff_beam)**i
        
        ## wavelength solution    
        dldp = np.zeros(NORDER)
        for i in range(NORDER):
            if 'DLDP_%s_%d' %(beam, i) in self.conf.keys():
                coeffs = self.conf['DLDP_%s_%d' %(beam, i)]
                dldp[i] = self.field_dependent(xi, yi, coeffs)
        
        ## dp is the arc length along the trace
        ## $\lambda = dldp_0 + dldp_1 dp + dldp_2 dp^2$ ...
        if self.conf['DYDX_ORDER_%s' %(beam)] == 0:   ## dy=0
            dp = dx-xoff_beam                      
        elif self.conf['DYDX_ORDER_%s' %(beam)] == 1: ## constant dy/dx
            dp = np.sqrt(1+dydx[1]**2)*(dx-xoff_beam)
        elif self.conf['DYDX_ORDER_%s' %(beam)] == 2: ## quadratic trace
            u0 = dydx[1]+2*dydx[2]*(0)
            dp0 = (u0*np.sqrt(1+u0**2)+np.arcsinh(u0))/(4*dydx[2])
            u = dydx[1]+2*dydx[2]*(dx-xoff_beam)
            dp = (u*np.sqrt(1+u**2)+np.arcsinh(u))/(4*dydx[2])-dp0
        else:
            ## high order shape, numerical integration along trace
            ## (this can be slow)
            xmin = np.minimum((dx-xoff_beam).min(), 0)
            xmax = np.maximum((dx-xoff_beam).max(), 0)
            xfull = np.arange(xmin, xmax)
            dyfull = 0
            for i in range(1, NORDER):
                dyfull += i*dydx[i]*(xfull-0.5)**(i-1)
            
            ## Integrate from 0 to dx / -dx
            dpfull = xfull*0.
            lt0 = xfull <= 0
            if lt0.sum() > 1:
                dpfull[lt0] = np.cumsum(np.sqrt(1+dyfull[lt0][::-1]**2))[::-1]
                dpfull[lt0] *= -1
            #
            gt0 = xfull >= 0
            if gt0.sum() > 0:
                dpfull[gt0] = np.cumsum(np.sqrt(1+dyfull[gt0]**2))
              
            dp = np.interp(dx-xoff_beam, xfull, dpfull)
        
        ## Evaluate dldp    
        lam = dp*0.
        for i in range(NORDER):
            lam += dldp[i]*dp**i
            
        return dy, lam
        
    def show_beams(self, beams=['E','D','C','B','A']):
        """
        Make a demo plot of the beams of a given configuration file
        """
        import matplotlib.pyplot as plt
        
        x0, x1 = 507, 507
        dx = np.arange(-800,1200)

        if 'WFC3.UV' in self.conf_file:
            x0, x1 = 2073, 250
            dx = np.arange(-1200,1200)
        if 'G800L' in self.conf_file:
            x0, x1 = 2124, 1024
            dx = np.arange(-1200,1200)
            
        s=200 # marker size
        fig = plt.figure(figsize=[10,3])
        plt.scatter(0,0,marker='s', s=s, color='black', edgecolor='0.8',
                    label='Direct')
        
        for beam in beams:
            if 'XOFF_%s' %(beam) not in self.conf.keys():
                continue
            
            xoff = self.field_dependent(x0, x1, self.conf['XOFF_%s' %(beam)])
            dy, lam = self.get_beam_trace(x0, x1, dx=dx, beam=beam)
            xlim = self.conf['BEAM%s' %(beam)]
            ok = (dx >= xlim[0]) & (dx <= xlim[1])
            plt.scatter(dx[ok]+xoff, dy[ok], c=lam[ok]/1.e4, marker='s', s=s,
                        alpha=0.5, edgecolor='None')
            plt.text(np.median(dx[ok]), np.median(dy[ok])+1, beam,
                     ha='center', va='center', fontsize=14)
            print 'Beam %s, lambda=(%.1f - %.1f)' %(beam, lam[ok].min(),
                                                    lam[ok].max())
            
        plt.grid()
        plt.xlabel(r'$\Delta x$')
        plt.ylabel(r'$\Delta y$')

        cb = plt.colorbar(pad=0.01, fraction=0.05)    
        cb.set_label(r'$\lambda\,(\mu\mathrm{m})$')
        plt.title(self.conf_file)
        plt.tight_layout()
        plt.savefig('%s.pdf' %(self.conf_file))
        
def get_extraction_range():
    """
    Compute good extraction range for 1st order spectra so that all 2D 
    spectra have same dimensions
    """
    grism_wlimit = {'G141':[1.05e4, 1.70e4, 22., 1.4e4], 'G102':[0.76e4, 1.17e4, 10., 1.05e4], 'G800L':[0.5e4, 1.05e4, 20., 0.75e4], 'GRS':[1.35e4, 1.95e4, 5., 1.65e4]}
    
    import mywfc3.grism
    
    grism = 'G141'
    
    conf = mywfc3.grism.aXeConf('%s.test27s.gbb.conf' %(grism))
    dx = np.arange(-50,400)
    xmi = 100
    xma = 0
    for ix in np.arange(0,1014,100):
        for iy in np.arange(0,1014,100):
            dy, wave = conf.get_beam_trace(x=ix, y=iy, dx=dx, beam='A')
            xminmax = np.interp(grism_wlimit[grism][0:2], wave, dx)
            print '%d %d  %f %f' %(ix, iy, xminmax[0], xminmax[1])
            xmi = np.minimum(xmi, xminmax[0])
            xma = np.maximum(xma, xminmax[1])
    
    grism = 'G800L'
    conf = mywfc3.grism.aXeConf('ACS.WFC.CHIP1.Cycle13.5.conf')
    conf = mywfc3.grism.aXeConf('ACS.WFC.CHIP2.Cycle13.5.conf')
    dx = np.arange(-50,400)
    xmi = 100
    xma = 0
    for ix in np.arange(0,4096,100):
        for iy in np.arange(0,2048,100):
            dy, wave = conf.get_beam_trace(x=ix, y=iy, dx=dx, beam='A')
            xminmax = np.interp(grism_wlimit[grism][0:2], wave, dx)
            print '%d %d  %f %f' %(ix, iy, xminmax[0], xminmax[1])
            xmi = np.minimum(xmi, xminmax[0])
            xma = np.maximum(xma, xminmax[1])
    
    ###
    first_order = {'G141':[28,178], 'G102':[45, 227], 'G800L':[60, 192]}
    
def compare_zeroth_orders():
    """
    Compare zeroth orders for different configuration files
    """
    import mywfc3.grism
    import unicorn
    
    root = 'FIGS-GN1-258'
    grow=2
    
    field = '-'.join(root.split('-')[:2])
    m = unicorn.reduce.GrismModel(root=root, growx=grow, growy=grow, direct='F105W', grism='G102')
    model_list = m.get_eazy_templates()
    #
    model = unicorn.reduce.process_GrismModel(root=root, grow_factor=grow, growx=grow, growy=grow, MAG_LIMIT=25, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, model_slope=0, model_list=model_list, direct='F105W', grism='G102', BEAMS=['A', 'B', 'C', 'D', 'E'], old_filenames=False)
        
    conf_list = {'G102.test27s.gbb.conf':'green', 'G102.test41.conf':'magenta'}
    conf_files = {}
    
    for f in conf_list.keys():
        conf_files[f] = mywfc3.grism.aXeConf('%s/CONF/%s' %(os.getenv('THREEDHST'), f))
    
    x_B = np.arange(-280,-220,0.1)
    l_B = 9100. # G102
    
    x_A = np.arange(41,248,0.1)
    l_A = [7900,1.135e4]
    

    ok = model.cat.mag < 22
    ix = np.arange(len(model.cat.mag))[ok]
    
    fp = open('/tmp/zeroth.reg','w')
    fp.write('image\n')
    
    for i in ix:
        xinter, yinter = model.cat.x_pix[i], model.cat.y_pix[i]
        xflt = (xinter-model.pad/2-model.ngrowx*model.growx)/model.growx
        yflt = (yinter-model.pad/2-model.ngrowy*model.growy)/model.growy
        
        for f in conf_list.keys():
            dy, lam = conf_files[f].get_beam_trace(x=xflt, y=yflt, dx=x_B, beam='B')
            xi = np.interp(l_B, lam, x_B)
            yi = np.interp(l_B, lam, dy)
            
            x0 = (xflt + xi)*model.growx + model.pad/2 + model.ngrowx*model.growx
            y0 = (yflt + yi)*model.growy + model.pad/2 + model.ngrowy*model.growy
            fp.write('circle(%.2f,%.2f,2) # color=%s width=2\n' %(x0, y0, conf_list[f]))
            #
            dy, lam = conf_files[f].get_beam_trace(x=xflt, y=yflt, dx=x_A, beam='A')
            xi = np.interp(l_A, lam, x_A)
            yi = np.interp(l_A, lam, dy)
            
            x0 = (xflt + xi)*model.growx + model.pad/2 + model.ngrowx*model.growx
            y0 = (yflt + yi)*model.growy + model.pad/2 + model.ngrowy*model.growy
            fp.write('line(%.2f,%.2f,%.2f,%.2f) # color=%s width=1\n' %(x0[0], y0[0], x0[1], y0[1], conf_list[f]))
            
            
    fp.close()
    
def go_extend():
    import mywfc3.grism
    
    for file in ['WFC3.IR.G102.2nd.sens.2.fits', 'WFC3.IR.G102.3rd.sens.2.fits', 'WFC3.IR.G141.2nd.sens.2.fits', 'WFC3.IR.G141.3rd.sens.2.fits']:
        mywfc3.grism.extend_sens(sens_file='/Users/brammer/3DHST/Spectra/Work/CONF/%s' %(file))
        
def extend_sens(sens_file='/Users/brammer/3DHST/Spectra/Work/CONF/WFC3.IR.G141.3rd.sens.2.fits'):
    """
    Some of the "v2" sensitivity files are cut off at the wavelength extremes for no apparent reason.  
    
    Here, extend them with multi-gaussian component fits 
    """
    
    import os
    import matplotlib.pyplot as plt
    
    from astroML.sum_of_norms import sum_of_norms, norm
    from threedhst import catIO
    

    #sens = catIO.Table('/Users/brammer/3DHST/Spectra/Work/CONF/WFC3.IR.G102.2nd.sens.2.fits')
    #sens = catIO.Table('/Users/brammer/3DHST/Spectra/Work/CONF/WFC3.IR.G102.3rd.sens.2.fits')
    if 'G102' in sens_file:
        xarr = np.arange(7100,1.19e4)
        n_gaussians, spacing = 20, 'linear'
    else:
        xarr = np.arange(9700,1.75e4)
        n_gaussians, spacing = 20, 'log'

    #sens = catIO.Table('/Users/brammer/3DHST/Spectra/Work/CONF/WFC3.IR.G141.2nd.sens.2.fits')
    sens = catIO.Table(sens_file)
    
    
    ok = sens['SENSITIVITY'] > 0
    
    x = sens['WAVELENGTH'][ok]
    y = sens['SENSITIVITY'][ok]
    ye = sens['ERROR'][ok]

    out = {}
    fig = plt.figure()
    
    for i, col in enumerate(['SENSITIVITY', 'ERROR']):
        ax = fig.add_subplot(211+i)
        y = sens[col][ok]
        
        w_best, rms, locs, widths = sum_of_norms(x, y, n_gaussians,
                                                 spacing=spacing,
                                                 full_output=True)

        norms = w_best * norm(xarr[:, None], locs, widths)

        out[col] = norms*1
        
        ax.plot(sens['WAVELENGTH'], sens[col], '-k', linewidth=3, label=col)
        ylim = ax.get_ylim()
        
        if col == 'SENSITIVITY':
            g = out[col].sum(1)
            keep = g > 1.e-4*g.max()
            ax.set_title(os.path.basename(sens.filename))
            
        ax.plot(xarr[keep], out[col][keep,:], ls='-', c='#FFAAAA')
        ax.plot(xarr[keep], out[col][keep,:].sum(1), '-r', linewidth=3, alpha=0.7)
        ax.set_ylim(-0.1 * ylim[1], 1.05*ylim[1])
        ax.legend(loc='upper right', fontsize=9)
    
    ax.set_xlabel('wavelength')
    fig.tight_layout(pad=1)
    
    t = catIO.table_base()
    t.add_column(catIO.Column(name='WAVELENGTH', data=xarr[keep]))
    for i, col in enumerate(['SENSITIVITY', 'ERROR']):
        t.add_column(catIO.Column(name=col, data=out[col][keep,:].sum(1)))
    
    outfile = sens.filename.replace('sens.2.fits', 'sens.X.fits')
    if os.path.exists(outfile):
        os.remove(outfile)
        
    t.write(outfile)
    print outfile
    plt.savefig(outfile.replace('.fits', '.png'))

    