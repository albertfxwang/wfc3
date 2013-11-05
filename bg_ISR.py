"""
Make figures for the ISR on backgrounds.
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
    
    fig = unicorn.plotting.plot_init(xs=12, aspect=0.45, left=0.1, top=0.05, bottom=0.1, hspace=0.0, wspace=0, use_tex=True, NO_GUI=True)

    ax, ax_limb = fig.add_subplot(256), fig.add_subplot(251)
    mywfc3.bg_ISR.orbit_tracks(FILTER='F105W', PATH='/Users/brammer/WFC3/Backgrounds/BroadBand/F105W', exposures=['ibp329iq', 'ibp329is'], axes=(ax, ax_limb), axlabels=[1,1,1])

    ax, ax_limb = fig.add_subplot(257), fig.add_subplot(252)
    mywfc3.bg_ISR.orbit_tracks(FILTER='F125W', PATH='/Users/brammer/WFC3/Backgrounds/BroadBand/Others', exposures=['ib5x19tm', 'ib5x19tq'], axes=(ax, ax_limb), axlabels=[1,0,0])
    
    ax, ax_limb = fig.add_subplot(258), fig.add_subplot(253)
    mywfc3.bg_ISR.orbit_tracks(FILTER = 'F160W',
    PATH = '/Users/brammer/WFC3/Backgrounds/BroadBand/Others',
    exposures = ['ib5x21ht', 'ib5x21hw'], axes=(ax, ax_limb), axlabels=[1,0,0])
    
    ax, ax_limb = fig.add_subplot(259), fig.add_subplot(254)
    mywfc3.bg_ISR.orbit_tracks(FILTER = 'G102',
    PATH = '/Users/brammer/WFC3/GrismPrograms/Stanford/RAW/',
    exposures = ['ibkn06dh', 'ibkn06dm', 'ibkn06dp', 'ibkn06dt'], axes=(ax, ax_limb), axlabels=[1,0,0])

    ax, ax_limb = fig.add_subplot(2,5,10), fig.add_subplot(255)
    # mywfc3.bg_ISR.orbit_tracks(FILTER = 'G141',
    # PATH = '/Users/brammer/WFC3/GrismPrograms/Koekemoer/RAW/',
    # exposures = ['ibl003a7', 'ibl003a9', 'ibl003aa', 'ibl003ab'],
    # axes=(ax, ax_limb), axlabels=[1,0,0])
    mywfc3.bg_ISR.orbit_tracks(FILTER = 'G141',
    PATH = '/Users/brammer/WFC3/Backgrounds/MultiAccum/',
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
    
    for i, exp in enumerate(exposures):
        dat = catIO.Readfile('%s/%sj_%s_orbit.dat' %(PATH, exp, FILTER), save_fits=False)
        raw = pyfits.open(gzfile('%s/%sq_raw.fits' %(PATH, exp)))
        bg_min = mywfc3.zodi.flt_zodi(raw.filename(), verbose=False)
        #
        spt = pyfits.getheader(gzfile('%s/%sq_spt.fits' %(PATH, exp)), 0)
        #### Start/stop times
        pstr = astropy.time.Time(spt['PSTRTIME'].replace('.',':'), format='yday', in_subfmt='date_hms', scale='utc')
        if i == 0:
           tstr = pstr
        #
        NSAMP = raw[0].header['NSAMP']
        times = np.zeros(NSAMP-1)
        for j in range(NSAMP-1):
            times[j] = raw['SCI', j+1].header['SAMPTIME']
        #    
        times = (times[::-1][1:] + (pstr-tstr).sec)/60.
        #
        ax.plot(times, dat.bg, color='black', linewidth=2)
        #
        day = dat.brightlimb > 0
        ax.plot(times[day], dat.bg[day], color='blue', alpha=0.4, linewidth=7, zorder=-10, label=['BrightLimb = 1','','',''][i])
        #
        l40 = dat.limbang > 40
        ax.plot(times[l40], dat.bg[l40], color='red', alpha=0.4, linewidth=4, zorder=-10, label=[r'LimbAng $>$ 40','','',''][i])
        #
        ax.plot(times, dat.bg*0.+bg_min, color='black', linestyle=':', linewidth=2, label=['Predicted zodi','','',''][i])
        ax.text(times[0]+0.2, 0.2+(N == 4)*((i % 2)*0.18-0.09), exp+'q', ha='left', va='bottom', fontsize=8)
        #
        ax_limb.plot(times, dat.limbang, color='black', linewidth=2)
        ax_limb.plot(times[day], dat.limbang[day], color='blue', alpha=0.4, linewidth=7, zorder=-10, label=['BrightLimb = 1','','',''][i])
        ax_limb.plot(times[l40], dat.limbang[l40], color='red', alpha=0.4, linewidth=4, zorder=-10, label=['BrightLimb = 1','','',''][i])
        ax_limb.plot(times, dat.termang, color='green', alpha=0.5, linewidth=2, label=['TermAng','','',''][i])
        
    ax.set_ylim(0,3.4)
    ax.set_xlim(-2,52)
    ax.set_xlabel(r'$\Delta t$ (minutes)')
    
    ax_limb.plot([-10,100], [20,20], linestyle='--', label='LimbAng = 20', color='black')
    ax_limb.set_ylim(0,100)
    ax_limb.set_xticklabels([])
    ax_limb.set_xlim(-2,52)
    ax_limb.set_title(FILTER)
    
    if axlabels[0] == 0:
        ax_limb.set_xticklabels([])

    if axlabels[1] == 0:
        ax_limb.set_yticklabels([])
        ax.set_yticklabels([])
    else:
        ax_limb.set_ylabel('Jitter: LimbAng')
        ax.set_ylabel('Background (electrons / s)')
        
    if axlabels[2]:
        ax.legend(loc='upper left', prop={'size':8})
        ax_limb.legend(loc='lower left', prop={'size':8})
        
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
    flat_file = os.path.join(os.getenv('iref'), '/Users/brammer/Research//HST/IREF/uc721143i_pfl.fits')
    
    direct = pyfits.open('../RAW/'+direct_flt)[1].data
    grism = []
    for file in grism_flt: grism.append(pyfits.open('../RAW/'+file)[1].data)
    flat = pyfits.open(flat_file)[1].data[5:-5,5:-5]
    
    #### Make flat-cleaned grism image
    #cleaned = 1./((grism[1]-grism[0])/flat)
    cleaned = flat/grism[1]
    excess = np.median(cleaned)
    cleaned -= excess
    
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
    
    red = unicorn.reduce.Interlace2D('red_00003.2D.fits')
    red.im['SCI'].data[:,:-1] = red.im['SCI'].data[:,1:]
    
    blue = unicorn.reduce.Interlace2D('blue_00002.2D.fits')
    
    er = np.loadtxt('red_excess.dat')
    eb = np.loadtxt('blue_excess.dat')
    
    xr, yr = red.optimal_extract(red.im['SCI'].data)
    xb, yb = blue.optimal_extract(blue.im['SCI'].data)
    
    xg = np.arange(5000,1.8e4,0.2)
    yg = np.exp(-(xg-1.0830e4)**2/2/2.**2)
    blue.compute_model(xg, yg)
    xb_model, yb_model = blue.optimal_extract(blue.model)

    red.compute_model(xg, yg)
    xr_model, yr_model = red.optimal_extract(red.model)
    
    fig = unicorn.plotting.plot_init(xs=8, aspect=0.6, left=0.01, top=0.01, bottom=0.01, right=0.01, hspace=0.0, wspace=0, use_tex=True, NO_GUI=True)
    
    ll = 0.06
    ax = fig.add_axes((ll, 0.1, 0.99-ll, 0.49))
    ax.plot(xb, eb/(yb+eb), color='blue', linewidth=2, alpha=0.8, label='G102, ibkn06dtq')
    ax.plot(xr, er/(yr+er), color='red', linewidth=2, alpha=0.8, label='G141, ibhj03xvq')
    ax.plot(1.083e4*np.ones(2), [0,1.5], color='green', linewidth=4, alpha=0.4, label=r'$\lambda = 10830\,$\AA$\,$ / convolved')
    ax.plot(xb_model, eb/(yb_model*4+eb), color='green', linewidth=2, alpha=0.5, zorder=10)
    ax.fill_between(xb_model, eb/(yb_model*4+eb), eb/(yb_model*4+eb)*0+1., color='green', linewidth=2, alpha=0.1, label=r'$\lambda = 10830\,$\AA, convolved', zorder=10, edgecolor='None')
    
    ax.set_xlim(0.783e4,1.383e4)
    xticks = np.array([0.8,0.9,1,1.1,1.2,1.3])*1.e4
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
        ax.set_xlim(scipy.polyval(a, [0.783e4, 1.383e4]))
        ax.set_xticks(scipy.polyval(a, xticks))
        ax.set_xticklabels([]); ax.set_yticklabels([]) 
        ax.set_yticks([0,2*NYSHOW])
        ax.set_ylabel(gname)
        
    unicorn.plotting.savefig(fig, 'extracted_blob_spectra.pdf')
    
    #### filter curves
    
    fig = unicorn.plotting.plot_init(xs=6, aspect=0.8, left=0.09, top=0.01, bottom=0.08, right=0.01, hspace=0.0, wspace=0., use_tex=True, NO_GUI=True)
    
    import pysynphot as S

    c = ['0.5','0.']
    ax = fig.add_subplot(311)
    for i, f in enumerate(['f098m','f125w']):
        bp = S.ObsBandpass('wfc3,ir,%s' %(f))
        ax.plot(bp.wave, bp.throughput, label=f.upper(), linewidth=2, alpha=0.8, color=c[i])
    #
    ax.set_xlim(0.783e4,1.38e4)
    ax.legend(loc='upper left')
    ax.set_xticklabels([])
    ax.plot(1.083e4*np.ones(2), [0,1.5], color='green', linewidth=4, alpha=0.4)
    ax.set_ylim(0,0.6)
    
    ax = fig.add_subplot(312)
    for i, f in enumerate(['f105w','f110w']):
        bp = S.ObsBandpass('wfc3,ir,%s' %(f))
        ax.plot(bp.wave, bp.throughput, label=f.upper(), linewidth=2, alpha=0.8, color=c[i])
        ax.fill_between(bp.wave, bp.throughput, bp.throughput*0, color='black', alpha=0.2)
        
    ax.set_xlim(0.783e4,1.38e4)
    ax.legend(loc='upper left')
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
    ax.legend(loc='upper left')
    ax.plot(1.083e4*np.ones(2), [0,1.5], color='green', linewidth=4, alpha=0.4)
    ax.set_ylim(0,0.6)
    ax.set_xlabel(r'$\lambda\,/\,$\AA')
    
    unicorn.plotting.savefig(fig, 'line_filter_transmissions.pdf')
    