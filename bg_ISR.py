"""
Make figures for the ISR on backgrounds.

ibsd02010 - nice track for figuring out scale height

"""
import matplotlib.pyplot as plt
import numpy as np
import pyfits

import astropy.time

import unicorn
from threedhst import catIO
import mywfc3.zodi
from mywfc3.utils import gzfile

def go_figure():
    import mywfc3.bg_ISR
    
    fig = unicorn.plotting.plot_init(xs=8, aspect=0.45, left=0.1, top=0.04, bottom=0.08, hspace=0.0, wspace=0, use_tex=True, NO_GUI=True)

    ax, ax_limb = fig.add_subplot(256), fig.add_subplot(251)
    mywfc3.bg_ISR.orbit_tracks(FILTER='F105W', PATH='/user/brammer/WFC3_Backgrounds/F105W/', exposures=['ibp329iq', 'ibp329is'], axes=(ax, ax_limb), axlabels=[1,1,0])

    ax, ax_limb = fig.add_subplot(257), fig.add_subplot(252)
    mywfc3.bg_ISR.orbit_tracks(FILTER='F125W', PATH='/user/brammer/WFC3_Backgrounds/F125W/', exposures=['ib5x19tm', 'ib5x19tq'], axes=(ax, ax_limb), axlabels=[1,0,1])
    
    ax, ax_limb = fig.add_subplot(258), fig.add_subplot(253)
    mywfc3.bg_ISR.orbit_tracks(FILTER = 'F160W',
    PATH = '/user/brammer/WFC3_Backgrounds/F160W/',
    exposures = ['ib5x21ht', 'ib5x21hw'], axes=(ax, ax_limb), axlabels=[1,0,0])
    
    ax, ax_limb = fig.add_subplot(259), fig.add_subplot(254)
    mywfc3.bg_ISR.orbit_tracks(FILTER = 'G102',
    PATH = '/user/brammer/WFC3_Backgrounds/GrismPrograms/',
    exposures = ['ibkn06dh', 'ibkn06dm', 'ibkn06dp', 'ibkn06dt'], axes=(ax, ax_limb), axlabels=[1,0,0])
    
    ax, ax_limb = fig.add_subplot(2,5,10), fig.add_subplot(255)
    # mywfc3.bg_ISR.orbit_tracks(FILTER = 'G141',
    # PATH = '/Users/brammer/WFC3/GrismPrograms/Koekemoer/RAW/',
    # exposures = ['ibl003a7', 'ibl003a9', 'ibl003aa', 'ibl003ab'],
    # axes=(ax, ax_limb), axlabels=[1,0,0])
    mywfc3.bg_ISR.orbit_tracks(FILTER = 'G141',
    PATH = '/user/brammer/WFC3_Backgrounds/GrismPrograms/',
    exposures = ['ibhj03xo', 'ibhj03xv'],
    axes=(ax, ax_limb), axlabels=[1,0,0])
    
    unicorn.plotting.savefig(fig, 'background_tracks.pdf')
    
def orbit_tracks(FILTER = 'F105W', PATH = '/Users/brammer/WFC3/Backgrounds/BroadBand/F105W', exposures = ['ibp329iq', 'ibp329is'], axes=None, axlabels=[1,1,1]):
    """
    Make a figure showing the orbit tracks / illuminated or not and the Term angle
    and compare to the observed background level (and predicted zodi).
    """
    
    ### F105W (all HUDF-DEEP-WFC3)
    # FILTER = 'F105W'
    # PATH = '/Users/brammer/WFC3/Backgrounds/BroadBand/F105W'
    # exposures = ['ibp329iq', 'ibp329is']
    
    # FILTER = 'F125W'
    # PATH = '/Users/brammer/WFC3/Backgrounds/BroadBand/Others'
    # exposures = ['ib5x19tm', 'ib5x19tq']
    # 
    # FILTER = 'F160W'
    # PATH = '/Users/brammer/WFC3/Backgrounds/BroadBand/Others'
    # exposures = ['ib5x21ht', 'ib5x21hw']
    
    # FILTER = 'G102'
    # PATH = '/Users/brammer/WFC3/GrismPrograms/Koekemoer/RAW/'
    # exposures = ['ibl003ac', 'ibl003ae', 'ibl003af', 'ibl003ag']
    # PATH = '/Users/brammer/WFC3/GrismPrograms/Stanford/RAW/'
    # exposures = ['ibkn06dh', 'ibkn06dm', 'ibkn06dp', 'ibkn06dt']

    # FILTER = 'G141'
    # PATH = '/Users/brammer/WFC3/GrismPrograms/Koekemoer/RAW/'
    # exposures = ['ibl003a7', 'ibl003a9', 'ibl003aa', 'ibl003ab']

    # PATH = '/Users/brammer/WFC3/GrismPrograms/CANDELS-SNe/RAW/'
    # exposures = ['ibfug1je', 'ibfug1jh']
    
    if axes is None:
        fig = unicorn.plotting.plot_init(xs=3, aspect=1/0.5, left=0.1, top=0.05, bottom=0.1, hspace=0.02, use_tex=True, NO_GUI=True)
        ax = fig.add_subplot(212)
        ax_limb = fig.add_subplot(211)
    else:
        ax, ax_limb = axes
    
    N = len(exposures)
    
    ax.plot([-10,100], [200,200], linestyle='--', label='LimbAng = 20', color='black') ## dummy for label
    
    for i, exp in enumerate(exposures):
        dat = catIO.Readfile('%s/%sj_%s_orbit.dat' %(PATH, exp, FILTER), save_fits=False)
        raw = pyfits.open(gzfile('%s/%sq_raw.fits' %(PATH, exp)))
        #bg_min = mywfc3.zodi.flt_zodi(raw.filename(), verbose=False)
        bg_min = dat.zodi
        #
        spt = pyfits.getheader(gzfile('%s/%sq_spt.fits' %(PATH, exp)), 0)
        #### Start/stop times
        pstr = astropy.time.Time(spt['PSTRTIME'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
        if i == 0:
           tstr = pstr
        
        NSAMP = raw[0].header['NSAMP']
        times = np.zeros(NSAMP-1)
        for j in range(NSAMP-1):
            times[j] = raw['SCI', j+1].header['SAMPTIME']
            
        times = (times[::-1][1:] + (pstr-tstr).sec)/60.
        
        ax.plot(times, dat.bg, color='black', linewidth=2)
        l40 = dat.limbang > 40
        ax.plot(times[l40], dat.bg[l40], color='red', alpha=0.4, linewidth=4, zorder=-10, label=r'LimbAng $>$ 40'*(i==0))
        
        shadow = dat.shadow == 1
        ax.plot(times[shadow], dat.bg[shadow], color='black', alpha=0.2, linewidth=7, zorder=-20, label=r'SHADOW'*(i==0))
        ax_limb.plot(times[shadow], dat.limbang[shadow], color='black', alpha=0.2, linewidth=7, zorder=-20, label=r'SHADOW'*(i==0))

        day = dat.brightlimb > 0
        ax.plot(times[day], dat.bg[day], color='blue', alpha=0.4, linewidth=7, zorder=-10, label='BrightLimb = 1'*(i==0))
        
        #
        ax.plot(times, dat.bg*0.+bg_min, color='black', linestyle=':', linewidth=2, label='Predicted zodi'*(i==0))
        ax.text(times[0]+0.2, 0.04+(N == 4)*((i % 2)*0.2-0.1), exp+'q', ha='left', va='bottom', fontsize=6)
        #
        ax_limb.plot(times, dat.limbang, color='black', linewidth=2)
        ax_limb.plot(times[day], dat.limbang[day], color='blue', alpha=0.4, linewidth=7, zorder=-10)
        ax_limb.plot(times[l40], dat.limbang[l40], color='red', alpha=0.4, linewidth=4, zorder=-10)
        #ax_limb.plot(times, dat.termang, color='green', alpha=0.5, linewidth=2, label=['TermAng','','',''][i])
        
    ax.set_ylim(-0.3,3.4)
    ax.set_xlim(-5,55)
    ax.set_xlabel(r'$\Delta\,t$ (minutes)')
    
    ax_limb.plot([-10,100], [20,20], linestyle='--', color='black')
    
    ax_limb.set_ylim(2,100)
    ax_limb.set_xticklabels([])
    ax_limb.set_xlim(-5,55)
    ax_limb.set_title(FILTER, size=10)
    
    if axlabels[0] == 0:
        ax_limb.set_xticklabels([])

    if axlabels[1] == 0:
        ax_limb.set_yticklabels([])
        ax.set_yticklabels([])
    else:
        ax_limb.set_ylabel('Jitter: LimbAng')
        ax.set_ylabel('Background (e- / s)')
        
    if axlabels[2]:
        ax.legend(loc='upper left', prop={'size':7})
        #ax_limb.legend(loc='lower left', prop={'size':8})
        
    if axes is None:
        unicorn.plotting.savefig(fig, 'track_%s.pdf' %(FILTER))
    
def extract_spectra():
    """
    Plots showing the line ID
    
    ### G102
    cp /Users/brammer/WFC3/GrismPrograms/Stanford/RAW/ibkn06*_[af]??.fits* /Users/brammer/WFC3/Backgrounds/FindLine/RAW
    
    ### G141
    cp /Users/brammer/WFC3/GrismPrograms/Koekemoer/RAW/ibl003*[af]??.fits* /Users/brammer/WFC3/Backgrounds/FindLine/RAW
    
    cp $THREEDHST/Incoming/ibhj03*[af]??.fits* /Users/brammer/WFC3/Backgrounds/FindLine/RAW
    """
    import os
    
    os.chdir('/Users/brammer/WFC3/Backgrounds/FindLine/Analysis')
    
    #### G102
    root, filt, gris = 'blue', 'F105W', 'G102'
    direct_flt = 'ibkn06d0q_flt.fits.gz'
    grism_flt = ['ibkn06dhq_flt.fits.gz', 'ibkn06dtq_flt.fits.gz']
    flat_file = os.path.join(os.getenv('iref'), '/Users/brammer/Research//HST/IREF//uc72113oi_pfl.fits')

    root, filt, gris = 'red', 'F140W', 'G141'
    direct_flt = 'ibhj03xtq_flt.fits.gz'
    grism_flt = ['ibhj03xoq_flt.fits.gz', 'ibhj03xvq_flt.fits.gz']

    grism_flt = ['ibhj03xoq_flt.fits', 'ibhj03xvq_flt.fits'] # EarthFlats/PREP

    flat_file = os.path.join(os.getenv('iref'), '/Users/brammer/Research//HST/IREF/uc721143i_pfl.fits')
    
    direct = pyfits.open('../RAW/'+direct_flt)[1].data
    grism = []
    for file in grism_flt: grism.append(pyfits.open('../RAW/'+file)[1].data)
    flat = pyfits.open(flat_file)[1].data[5:-5,5:-5]
    
    #### Make flat-cleaned grism image
    #cleaned = 1./((grism[1]-grism[0])/flat)
    cleaned = flat/grism[1]
    excess = np.median(cleaned)
    cleaned = cleaned/excess-1
    
    np.savetxt(root+'_excess.dat', [excess])
    
    im = pyfits.open('../RAW/'+grism_flt[1])
    im[1].data = cleaned
    im.writeto(grism_flt[1].split('.gz')[0], clobber=True)
    
    asn = threedhst.utils.ASNFile('../RAW/ibkn06030_asn.fits')
    asn.exposures = [grism_flt[1].split('_flt')[0]]
    asn.product = root+'-'+gris
    asn.write(root+'-'+gris+'_asn.fits')
    
    #### Make object image with the inverse blobs
    im = pyfits.open('../RAW/'+direct_flt)
    blobs = [(230, 418), (467, 377)]
    source = flat*0.
    for blob in blobs:
        yi, xi = np.indices(flat.shape)
        r = np.sqrt((xi-blob[0])**2+(yi-blob[1])**2)
        source[r < 20] += 1./flat[r < 20]-1
    
    im[1].data = source*100
    im.writeto(direct_flt.split('.gz')[0], clobber=True)
    
    asn.exposures = [direct_flt.split('_flt')[0]]
    asn.product = root+'-'+filt
    asn.write(root+'-'+filt+'_asn.fits')
        
    unicorn.reduce.interlace_combine(root+'-'+filt, growx=1, growy=1)
    unicorn.reduce.interlace_combine(root+'-'+gris, growx=1, growy=1)
    
    model = unicorn.reduce.GrismModel(root=root, grow_factor=1, growx=1, growy=1, MAG_LIMIT=20, use_segm=False, grism=gris, direct=filt)
    for obj in model.objects:
        model.twod_spectrum(obj, miny=-100, refine=True)
        model.show_2d(savePNG=True)
        
def find_line():
    import unicorn
    import scipy.ndimage as nd
    import scipy
    import os
    
    os.chdir('/Users/brammer/WFC3/Backgrounds/FindLine/Analysis/')
    
    #red = unicorn.reduce.Interlace2D('red_00002.2D.fits') # EarthFlats
    red = unicorn.reduce.Interlace2D('red_00003.2D.fits')

    red.im['SCI'].data[:,:-1] = red.im['SCI'].data[:,1:]
    
    blue = unicorn.reduce.Interlace2D('blue_00002.2D.fits')
    
    er = np.loadtxt('red_excess.dat')
    eb = np.loadtxt('blue_excess.dat')
    
    xr, yr = red.optimal_extract(red.im['SCI'].data)
    xb, yb = blue.optimal_extract(blue.im['SCI'].data)
    
    l0 = 1.0830e4
    xg = np.arange(5000,1.8e4,0.2)
    yg = np.exp(-(xg-l0)**2/2/2.**2)
    
    red.compute_model(xg, yg)
    xr_model, yr_model = red.optimal_extract(red.model)

    blue.compute_model(xg, yg)
    xb_model, yb_model = blue.optimal_extract(blue.model)
    
    #fig = unicorn.plotting.plot_init(xs=8, aspect=0.6, left=0.01, top=0.01, bottom=0.01, right=0.01, hspace=0.0, wspace=0, use_tex=True, NO_GUI=True)
    
    fig = unicorn.plotting.plot_init(xs=4, aspect=1, left=0.09, top=0.01, bottom=0.08, right=0.01, hspace=0.03, wspace=0., use_tex=True, NO_GUI=True)
    
    ll = 0.11
    ax = fig.add_axes((ll, 0.1, 0.99-ll, 0.49))
    ax.plot(xb, eb/(yb+eb), color='blue', linewidth=2, alpha=0.8, label='G102, ibkn06dtq')
    ax.plot(xr, er/(yr+er), color='red', linewidth=2, alpha=0.8, label='G141, ibhj03xvq')
    ax.plot(1.083e4*np.ones(2), [0,1.5], color='green', linewidth=4, alpha=0.4, label=r'$\lambda = 10830\,$\AA$\,$ / convolved')
    ax.plot(xb_model, eb/(yb_model*4+eb), color='green', linewidth=2, alpha=0.5, zorder=10)
    ax.fill_between(xb_model, eb/(yb_model*4+eb), eb/(yb_model*4+eb)*0+1., color='green', linewidth=2, alpha=0.1, label=r'$\lambda = 10830\,$\AA, convolved', zorder=10, edgecolor='None')
    
    xlim = [0.783e4, 1.383e4]
    xticks = np.array([0.8,0.9,1,1.1,1.2,1.3])*1.e4

    xlim = [0.933e4, 1.233e4]
    xticks = np.array([1., 1.1, 1.2])*1.e4
    
    ax.set_xlim(xlim)
    ax.set_xticks(xticks)
    ax.set_ylim(0,1.5)
    ax.set_xlabel(r'$\lambda\,/\,$\AA')
    ax.set_ylabel(r'Blob spectrum')
    ax.legend(loc='lower left', prop={'size':8})
    
    #### Show spectra
    NYSHOW = 20
    dy = (0.49-0.1)/2.
    bg = [eb, er]
    
    for i, spec, bg, wave, gname in zip([0,1], [red, blue], [er, eb], [xr, xb], ['G141', 'G102']):
        NY, NX = spec.im['SCI'].data.shape
        ax = fig.add_axes((ll, 0.1+0.49+dy*i, 0.99-ll, dy))
        spec2D = spec.im['SCI'].data
        ax.imshow(bg/(spec2D[NY/2-NYSHOW:NY/2+NYSHOW,:]+bg), aspect='auto', interpolation='nearest', vmin=0.8, vmax=1.1)
        xarr = np.arange(NX)
        a = scipy.polyfit(wave, xarr,1)
        ax.set_xlim(scipy.polyval(a, xlim))
        ax.set_xticks(scipy.polyval(a, xticks))
        ax.set_xticklabels([]); ax.set_yticklabels([]) 
        ax.set_yticks([0,2*NYSHOW])
        ax.set_ylabel(gname)
        
    unicorn.plotting.savefig(fig, 'extracted_blob_spectra.pdf')
    
    
def throughput_curves():
    #### filter curves
    fig = unicorn.plotting.plot_init(xs=4, aspect=1, left=0.09, top=0.01, bottom=0.08, right=0.01, hspace=0.03, wspace=0., use_tex=True, NO_GUI=True)
    
    import pysynphot as S

    c = ['0.5','0.']
    ax = fig.add_subplot(311)
    for i, f in enumerate(['f098m','f125w']):
        bp = S.ObsBandpass('wfc3,ir,%s' %(f))
        ax.plot(bp.wave, bp.throughput, label=f.upper(), linewidth=2, alpha=0.8, color=c[i])
    #
    ax.set_xlim(0.783e4,1.38e4)
    ax.legend(loc='upper left', prop={'size':9})
    ax.set_xticklabels([])
    ax.plot(1.083e4*np.ones(2), [0,1.5], color='green', linewidth=4, alpha=0.4)
    ax.set_ylim(0,0.6)
    
    ax = fig.add_subplot(312)
    for i, f in enumerate(['f105w','f110w']):
        bp = S.ObsBandpass('wfc3,ir,%s' %(f))
        ax.plot(bp.wave, bp.throughput, label=f.upper(), linewidth=2, alpha=0.8, color=c[i])
        ax.fill_between(bp.wave, bp.throughput, bp.throughput*0, color='black', alpha=0.2)
        
    ax.set_xlim(0.783e4,1.38e4)
    ax.legend(loc='upper left', prop={'size':9})
    ax.set_xticklabels([])
    ax.plot(1.083e4*np.ones(2), [0,1.5], color='green', linewidth=4, alpha=0.4)
    ax.set_ylim(0,0.6)
    ax.set_ylabel('Synphot throughput')
    
    ax = fig.add_subplot(313)
    for i, f in enumerate(['g102','g141']):
        bp = S.ObsBandpass('wfc3,ir,%s' %(f))
        ax.plot(bp.wave, bp.throughput, label=f.upper(), linewidth=2, alpha=0.8, color=c[i])
        ax.fill_between(bp.wave, bp.throughput, bp.throughput*0, color='black', alpha=0.2)
        
    ax.set_xlim(0.783e4,1.38e4)
    ax.legend(loc='upper left', prop={'size':9})
    ax.plot(1.083e4*np.ones(2), [0,1.5], color='green', linewidth=4, alpha=0.4)
    ax.set_ylim(0,0.6)
    ax.set_xlabel(r'$\lambda\,/\,$\AA')
    
    unicorn.plotting.savefig(fig, 'line_filter_transmissions.pdf')
    
def blob_demo():
    """
    Show the flat-field and spectra with and without the blob
    """
    import pyfits
    CONF = '/Users/brammer/3DHST/Spectra/Work/CONF/'
    no_blob = pyfits.open(CONF + 'sky.G141.set003.fits')[0].data
    yes_blob = pyfits.open(CONF + 'sky.G141.set025.fits')[0].data
    flat = pyfits.open('/Users/brammer/Research/HST/IREF/flat_3DHST_F140W_t1_v0.1.fits')[1].data[5:-5,5:-5]
    
    xy = np.array([[204,377], [431,457]])
    
    NX, NY = np.diff(xy, axis=0)[0]
    
    la, ra, ta, ba = 0.055, 0.005, 0.16, 0.29

    aspect = (NY*(1+ta+ba))/((3.*NX)*(1+la+ra))
    
    fig = unicorn.plotting.plot_init(xs=8, aspect=aspect, left=0, right=0, top=0, bottom=0, hspace=0.0, wspace=0, use_tex=True, NO_GUI=True)
    
    images = [flat, no_blob, yes_blob]
    labels = ['Flat F140W', 'G141/Flat, low background', 'G141/Flat, high background']
    tag = ['a)', 'b)', 'c)']
    dx = (1-la-ra)/3.    
    for i in range(3):
        ax = fig.add_axes((la+dx*i, ba, dx, 1-ba-ta))
        med = np.median(images[i][xy[1][0]:xy[1][1],xy[0][0]:xy[0][1]])
        print med
        ax.imshow(images[i]/med, vmin=0.92, vmax=1.04, aspect='auto')
        ax.set_xlim(xy[:,0])
        ax.set_ylim(xy[:,1])
        ax.set_yticks([400,450])
        if i > 0:
            ax.set_yticklabels([])
            ax.set_xlabel(r'$x\ (\lambda\rightarrow)$')
        else:
            ax.set_ylabel(r'$y$')
            ax.set_xlabel(r'$x$')
        #
        ax.set_title(labels[i], size=10)
        ax.text(0.0165, 0.94, tag[i], ha='left', va='top', backgroundcolor='white', size=8, transform=ax.transAxes)
        #
        # if i == 0:
        #     dxa, dya, dd = 15, -5, 20 
        #     thet = np.arctan2(dya, dxa)
        #     ax.arrow(230+dxa+dd*np.cos(thet), 417+dya+dd*np.sin(thet), -dd*np.cos(thet), -dd*np.sin(thet), color='white', linewidth=1.2, head_width=dd*0.2, head_length=dd*0.2, overhang=0.5, length_includes_head=False)
        #     ax.text(230+dxa+dd*np.cos(thet)+5, 417+dya+dd*np.sin(thet), 'Blob', color='black', ha='left', va='center')
        
    #
    unicorn.plotting.savefig(fig, 'blob_demo.pdf')
    
def excess_statistics():
    """
    Show cumulative distribution of (read) excess backgrounds
    """

    master = {}
    master['F105W'] = catIO.Readfile('/user/brammer/WFC3_Backgrounds/F105W/master.dat')
    master['G141'] = catIO.Readfile('/user/brammer/WFC3_Backgrounds/GrismPrograms/master_G141.dat')
    master['G102'] = catIO.Readfile('/user/brammer/WFC3_Backgrounds/GrismPrograms/master_G102.dat')
    
    la, ra, ta, ba = 0.06, 0.06, 0.07, 0.13

    NX, NY = 1, 1
    aspect = (NY*(1+ta+ba))/((3.*NX)*(1+la+ra))
    dx = (1-la-ra)/3.    
    
    fig = unicorn.plotting.plot_init(xs=8, aspect=aspect, left=0, right=0, top=0, bottom=0, hspace=0.0, wspace=0, use_tex=True, NO_GUI=True)

    for i in range(3):
        bgf = master[master.keys()[i]]
        ax = fig.add_axes((la+dx*i, ba, dx, 1-ba-ta))
        #
        bg_ratio, xr = bgf.bg-bgf.zodi, (-0.5,6)
        bg_ratio, xr = bgf.bg/bgf.zodi, (0.5,6.8)
        #
        yh1, xh1, nn = ax.hist(bg_ratio[bgf.shadow == 1], range=xr, bins=100, alpha=0.2, color='black', histtype='stepfilled', log=True)
        yhc1 = np.cumsum(yh1[::-1])
        #
        yh0, xh0, nn = ax.hist(bg_ratio[bgf.shadow == 0], range=xr, bins=100, alpha=0.2, color='red', histtype='stepfilled', log=True)
        yhc0 = np.cumsum(yh0[::-1])
        #
        axr = ax.twinx()
        #
        if i == 0:
            ax.plot([100,110],[1.e4,1.e4], linewidth=8, color='black', alpha=0.2, label='SHADOW = True')
            ax.plot([100,110],[1.e4,1.e4], linewidth=8, color='red', alpha=0.2, label='SHADOW = False')
            ax.plot([100,110],[1.e4,1.e4], linewidth=1.2, color='red', alpha=1, label=r'$f(>X)$')
            ax.legend(loc='right', prop={'size':9}, scatterpoints=1, frameon=True)
            
        #
        yi = np.interp(2, xh0[1:], yhc0[::-1]*1./bgf.N)
        axr.scatter([2],[yi], color='red', marker='o', label=r'$f(X>2)$ = %.1f' %(yi*100) + '\%')
        axr.plot(xh1[1:][::-1], yhc1*1./bgf.N, color='black', linewidth=1.2)
        axr.plot(xh0[1:][::-1], yhc0*1./bgf.N, color='red', linewidth=1.2)
        #
        ### efficiency = 1/X
        eff = 1/(np.clip(bgf.bg/bgf.zodi,1,10))
        axr.scatter([100,110],[-1,-1], linewidth=1, color='white', alpha=0.2, label='eff. = %d' %(np.sum(eff)/bgf.N*100) + '\%')
        
        axr.legend(loc='upper right', prop={'size':9}, scatterpoints=1, frameon=False)
        ax.set_xlim(xr)
        #ax.set_ylim(5,yh1.max()*1.2)
        ax.set_ylim(5,3000)
        axr.set_ylim(0,0.65)
        if i > 0:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel(r'$N_\mathrm{read}$')
        if i < 2:
            axr.set_yticklabels([])
        else:
            axr.set_ylabel(r'$f(>X)$')
        
        ax.set_title(master.keys()[i], size=10)
        if i == 1:
            ax.set_xlabel(r'$X$ = background observed / predicted zodi')
    #
    unicorn.plotting.savefig(fig, '/tmp/excess_statistics.pdf')
    