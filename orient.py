"""
Compute POS-TARG offsets for maximum overlap of 100% spectra
"""

import numpy as np
import matplotlib.pyplot as plt

from shapely.geometry import Polygon, Point, MultiPolygon, LineString
from descartes import PolygonPatch

import unicorn
import threedhst

def overlap(theta=90, grism='G102', have_mosaic=True, plot=True, aper='GRISM1024', f_spec=1, ax=None, get_postargs=False, recenter_target=False):
            
    refpix_aper = {'GRISM1024':[497, 562], 'IR':[562, 562], 'IR-FIX': [512,512]}
    disp_box_aper = {'G102': [53, 210], 'G141': [36, 168]}

    #### Pixel scale
    ps = {'x':0.1355, 'y':0.1211} 
    
    refpix = np.array(refpix_aper[aper])
    disp_box = np.array(disp_box_aper[grism])
    
    spec_width = np.diff(disp_box)[0]
    #disp_box[1] -= spec_width*(1-f_spec) + disp_box[0]
    right = 1014 - disp_box[0] - spec_width*(f_spec)
    
    #### Assume have direct imaging at 
    if have_mosaic:
        left = (-1*disp_box[0]) - spec_width*(1-f_spec)
    else:
        left = 0 #- spec_width*(1-f_spec)
    
    #### Full box: full WFC3/IR
    full_box_x = np.array([0, 1014, 1014, 0, 0])*ps['x']
    full_box_y = np.array([0, 0, 1014, 1014, 0])*ps['y']
    
    #sub_box_x = np.array([disp_box[0], 1014-disp_box[1], 1014-disp_box[1], disp_box[0], disp_box[0]])*ps['x']
    sub_box_x = np.array([left, right, right, left, left])*ps['x']
    sub_box_y = full_box_y*1.
    
    #### Center of rotation of sub box
    sub_center_x = ((right-left)/2.+left)*ps['x']
    sub_center_y = 1014./2*ps['y']
    
    #### Center of rotation for aperture
    aper_center_x = refpix[0]*ps['x']
    aper_center_y = refpix[1]*ps['y']
    
    #### Compute POS-TARGs
    if recenter_target:
        recenter_x = aper_center_x - sub_center_x
        recenter_y = aper_center_y - sub_center_y        
        
        print 'Recenter POSTARGs: %.3f %.3f' %(-recenter_x, -recenter_y)
        
        full_box_x += recenter_x
        full_box_y += recenter_y

        sub_box_x += recenter_x
        sub_box_y += recenter_y
    
        
    full_poly = Polygon(np.array([full_box_x, full_box_y]).T)
    full_patch = PolygonPatch(full_poly, fc='None', ec='black', alpha=0.5, zorder=2)
    
    sub_poly = Polygon(np.array([sub_box_x, sub_box_y]).T)
    sub_patch = PolygonPatch(sub_poly, fc='black', ec='black', alpha=0.5, zorder=2)
    
    #### Rotate about aperture center
    rot_full_aper_x, rot_full_aper_y = threedhst.utils.xyrot(full_box_x, full_box_y, theta, x0=aper_center_x, y0=aper_center_y, radians=False, ccw=False)
    rot_sub_aper_x, rot_sub_aper_y = threedhst.utils.xyrot(sub_box_x, sub_box_y, theta, x0=aper_center_x, y0=aper_center_y, radians=False, ccw=False)
    
    rot_full_aper_poly = Polygon(np.array([rot_full_aper_x, rot_full_aper_y]).T)
    rot_full_aper_patch = PolygonPatch(rot_full_aper_poly, fc='None', ec='red', alpha=0.5, zorder=2)
    rot_sub_aper_poly = Polygon(np.array([rot_sub_aper_x, rot_sub_aper_y]).T)
    rot_sub_aper_patch = PolygonPatch(rot_sub_aper_poly, fc='red', ec='red', alpha=0.5, zorder=2)
    
    aper_isect = sub_poly.intersection(rot_sub_aper_poly)
    
    #### Rotate about coverage center
    rot_full_cov_x, rot_full_cov_y = threedhst.utils.xyrot(full_box_x, full_box_y, theta, x0=sub_center_x, y0=sub_center_y, radians=False, ccw=False)
    rot_sub_cov_x, rot_sub_cov_y = threedhst.utils.xyrot(sub_box_x, sub_box_y, theta, x0=sub_center_x, y0=sub_center_y, radians=False, ccw=False)
    
    rot_full_cov_poly = Polygon(np.array([rot_full_cov_x, rot_full_cov_y]).T)
    rot_full_cov_patch = PolygonPatch(rot_full_cov_poly, fc='None', ec='blue', alpha=0.5, zorder=2)
    rot_sub_cov_poly = Polygon(np.array([rot_sub_cov_x, rot_sub_cov_y]).T)
    rot_sub_cov_patch = PolygonPatch(rot_sub_cov_poly, fc='blue', ec='blue', alpha=0.5, zorder=2)

    cov_isect = sub_poly.intersection(rot_sub_cov_poly)
    
    #### POS TARGs
    pt_x_rot = (rot_full_cov_x - rot_full_aper_x)#[0]
    pt_y_rot = (rot_full_cov_y - rot_full_aper_y)#[0]
    
    pt_x, pt_y = threedhst.utils.xyrot(pt_x_rot, pt_y_rot, -theta, x0=0, y0=0, ccw=False, radians=False)
    
    if get_postargs:
        return pt_x[0], pt_y[0]
        
    fig = None
    if plot:
        if ax is None:
            fig = unicorn.plotting.plot_init(xs=6, aspect=1., square=True, left=0.1, bottom=0.1, right=0.01, top=0.01)
        #ax = fig.add_subplot(111)
            ax = fig.add_axes((0.1,0.1,0.89, 0.89))
        
        ax.scatter(-200,-200)
        ax.set_xlim(-150,200)
        ax.set_ylim(-150,200)
        ax.add_patch(full_patch)
        ax.add_patch(sub_patch)

        ax.add_patch(rot_full_aper_patch)
        ax.add_patch(rot_sub_aper_patch)

        ax.add_patch(rot_full_cov_patch)
        ax.add_patch(rot_sub_cov_patch)
        
        ax.scatter(aper_center_x, aper_center_y, marker='x', color='red', label='APER center (%5.3f)' %(aper_isect.area/full_poly.area), zorder=10)
        ax.scatter(sub_center_x, sub_center_y, marker='o', color='blue', label='Coverage center (%5.3f)' %(cov_isect.area/full_poly.area), zorder=10)
        
        ax.legend(loc='lower right', title=r'$\theta$ = %.2f, $f$=%.2f' %(theta, f_spec), scatterpoints=1)

        #plt.draw()
        
    return rot_full_aper_poly, rot_sub_aper_poly, aper_isect.area/full_poly.area, rot_full_cov_poly, rot_sub_cov_poly, cov_isect.area/full_poly.area, fig
    
def show_all():
    import mywfc3.orient
    
    thetas = np.arange(-180, 181)
    area_aper = thetas*0.
    area_cov = thetas*0.
    
    f_spec = 0.75
    for i in thetas:
        result = mywfc3.orient.overlap(theta=thetas[i], grism='G102', have_mosaic=False, plot=False, aper='GRISM1024', f_spec=f_spec)
        area_aper[i], area_cov[i] = result[2], result[5] 
        
    plt.plot(thetas, area_aper, label='APER', color='red')
    plt.plot(thetas, area_cov, label='Recenter', color='blue')
    
    #### FIGs
    ORIENTS = [347, 279, 217, 208, 286] # A
    
    so = np.argsort(ORIENTS)
    ix = np.cast[int](ORIENTS-np.median(ORIENTS))[so] + 180
    
    plt.scatter(thetas[ix], area_aper[ix], color='red')
    plt.scatter(thetas[ix], area_cov[ix], color='blue')
    for i in so:
        ii = np.cast[int](ORIENTS-np.median(ORIENTS))[i] + 180
        plt.text(thetas[ii], area_aper[ii]-0.02, ORIENTS[i], fontsize=8, ha='center', va='top')
        
    plt.legend(loc='lower left', title=r'$f$ > %.2f' %(f_spec))
    plt.xlabel(r'$\theta = \Delta$ ORIENT (deg)')
    plt.ylabel('Area with > %2d%s trace coverage' %(f_spec*100,'%'))
    plt.grid()
    plt.savefig('orient_area_fraction.pdf'); plt.close()
    
    result = mywfc3.orient.overlap(theta=68, grism='G102', have_mosaic=True, plot=True, aper='GRISM1024', f_spec=f_spec)
    unicorn.plotting.savefig(result[-1], 'orient_overlap+68.pdf'); plt.close()
    
    result = mywfc3.orient.overlap(theta=-71, grism='G102', have_mosaic=True, plot=True, aper='GRISM1024', f_spec=f_spec)
    unicorn.plotting.savefig(result[-1], 'orient_overlap-71.pdf'); plt.close()
    
    postarg_x = thetas*0.
    postarg_y = thetas*0.
    for i in thetas:
        postarg_x[i], postarg_y[i] = mywfc3.orient.overlap(theta=thetas[i], grism='G102', have_mosaic=True, plot=False, aper='GRISM1024', f_spec=f_spec, get_postargs=True)
    
    #
    plt.plot(thetas, postarg_x, label='POSTARG X', color='black')
    plt.plot(thetas, postarg_y, label='POSTARG Y', color='green')
    
    plt.scatter(thetas[ix], postarg_x[ix], color='black')
    plt.scatter(thetas[ix], postarg_y[ix], color='green')
    for i in so:
        ii = np.cast[int](ORIENTS-np.median(ORIENTS))[i] + 180
        plt.text(thetas[ii], postarg_x[ii]-1, ORIENTS[i], fontsize=8, ha='center', va='top')
        
    plt.grid()
    plt.legend(loc='lower left', scatterpoints=1)
    plt.xlabel(r'$\theta = \Delta$ ORIENT (deg)')
    plt.ylabel('POSTARG (arcsec)')
    plt.savefig('orient_align_postargs.pdf'); plt.close()
    
    #### Exposure maps
    results = []
    for theta in np.cast[int](ORIENTS-np.median(ORIENTS))[so]:
        result = mywfc3.orient.overlap(theta=theta, grism='G102', have_mosaic=True, plot=False, aper='GRISM1024', f_spec=f_spec, recenter_target=False)
        results.append(result)
    #
    # results = []
    # for theta in np.arange(0,361,10):
    #     result = mywfc3.orient.overlap(theta=theta, grism='G102', have_mosaic=False, plot=False, aper='GRISM1024', f_spec=f_spec, recenter_target=False)
    #     results.append(result)
    
    alpha = 1./len(results)*1.5
    ### No rotation
    j=0
    
    fig = unicorn.plotting.plot_init(xs=10, aspect=0.5, square=True, left=0.01, bottom=0.01, right=0.01, top=0.01)

    if j == 0:
        ax = fig.add_axes((0.01,0.01,0.48, 0.98))
        rot = 'APER'
        full_depth = results[0][1].intersection(results[0][1])
    else:
        ax = fig.add_axes((0.51,0.01,0.48, 0.98))
        rot = 'Recenter'
        full_depth = results[0][4].intersection(results[0][4])

    ax.scatter(-200,-200); ax.set_xlim(-80,200); ax.set_ylim(-80,200)
    ax.set_xticklabels([]); ax.set_yticklabels([])
        
    for result in results:
        if rot == 'APER':
            full_poly, sub_poly = result[0], result[1]
            full_patch = PolygonPatch(full_poly, fc='None', ec='red', alpha=alpha, zorder=2)
            sub_patch = PolygonPatch(sub_poly, fc='red', ec='red', alpha=alpha, zorder=2)
        else:
            full_poly, sub_poly = result[3], result[4]
            full_patch = PolygonPatch(full_poly, fc='None', ec='blue', alpha=alpha, zorder=2)
            sub_patch = PolygonPatch(sub_poly, fc='blue', ec='blue', alpha=alpha, zorder=2)
        #   
        ax.add_patch(full_patch)
        ax.add_patch(sub_patch)
        #
        full_depth = full_depth.intersection(sub_poly)
    
    full_depth_patch = PolygonPatch(full_depth, fc='None', ec='black', alpha=0.8, zorder=2)
    ax.add_patch(full_depth_patch)
    ax.plot([-900],[-900], color='black', label=r'Full depth: %.2f arcmin$^2$' %(full_depth.area/3600.))
    ax.set_xlim(-80,200); ax.set_ylim(-80,200)
    
    ax.legend(loc='lower left', title='Rotation: %s' %(rot))
    
    unicorn.plotting.savefig(fig, 'orient_full_depth.pdf')
    
    result = mywfc3.orient.overlap(theta=0, grism='G102', have_mosaic=True, plot=False, aper='GRISM1024', f_spec=f_spec, recenter_target=True)
    
def ISR_figures_demo():
    
    import mywfc3.orient
    #### Pixel scale
    ps = {'x':0.1355, 'y':0.1211} 
    #### dispersion offsets
    disp_box = {'G102': [53, 210], 'G141': [36, 168]}
    #### Apertures
    refpix_aper = {'GRISM1024':[497, 562], 'IR':[562, 562], 'IR-FIX': [512,512]}
    
    grism = 'G102'
    aper = 'GRISM1024'
    
    #### Demo
    fig = unicorn.plotting.plot_init(xs=11, aspect=1./4, square=True, left=0.01, bottom=0.01, right=0.01, top=0.01, use_tex=True, NO_GUI=True)
    
    theta_pairs = [[-20], [50], [160], [-20,50,160]]
    
    dx = 0.99/len(theta_pairs)
    
    for i in range(len(theta_pairs)):
        
        ax = fig.add_axes((0.005+i*(dx+0*0.005),0.01, dx, 0.97))    
        
        if i < 3:
            ax.text(0.5, 0.05, r'PA($y$) = %d' %(-theta_pairs[i][0]), transform=ax.transAxes, ha='center', va='bottom')
        else:
            ax.text(0.5, 0.05, r'PA($y$) = %d $\times$ %d $\times$ %d' %(-theta_pairs[i][0], -theta_pairs[i][1], -theta_pairs[i][2]), transform=ax.transAxes, ha='center', va='bottom')
            
        ax.scatter(refpix_aper[aper][0]*ps['x'],refpix_aper[aper][1]*ps['y'], color='purple', marker='x', s=60, label=aper)
        plt_range = np.array([-72, 180])*0.95
        ax.set_xlim(plt_range+10); ax.set_ylim(plt_range+10)
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.set_xticks(ax.get_xlim()); ax.set_yticks(ax.get_ylim())
        
        ax.plot([-1000,-1100],[-1110,-1100], color='green', alpha=0.4, linewidth=2, label=grism)
        
        if i == 0:
            ax.legend(scatterpoints=1, loc='upper left', fontsize=10)
            
        #theta = 70
        pair = theta_pairs[i]
        results = []
        for theta in pair:
            result = mywfc3.orient.overlap(theta=theta, grism=grism, have_mosaic=True, plot=False, aper='GRISM1024', f_spec=0.5, recenter_target=False)
            results.append(result)
        
        overlap_region = results[0][1].intersection(results[0][1])
        for result in results:
            overlap_region = overlap_region.intersection(result[1])
            
        for j in range(len(pair)):
            full_poly, sub_poly = results[j][0], results[j][1]            
            full_patch = PolygonPatch(full_poly, fc='None', ec='black', alpha=0.9, zorder=-2)
            sub_patch = PolygonPatch(sub_poly, fc='black', ec='black', alpha=0.2, zorder=-2)
            
            ax.add_patch(full_patch)
            ax.add_patch(sub_patch)
                        
            ### Death star
            xds, yds = threedhst.utils.xyrot(np.array([357,357])*ps['x'], np.array([55,55])*ps['y'], pair[j], x0=refpix_aper[aper][0]*ps['x'], y0=refpix_aper[aper][1]*ps['y'])
            ds = Point(xds[0], yds[0]).buffer(20*ps['x'])
            ds_patch = PolygonPatch(ds, fc='None', ec='black', alpha=0.9, zorder=-2)
            ax.add_patch(ds_patch)
            
            if i == 3:
                isect_patch = PolygonPatch(overlap_region, fc='green', ec='green', alpha=0.2, zorder=1)
                ax.add_patch(isect_patch)
                
            else:
                for xi in np.arange(-15, 180, 33):
                    for yi in np.arange(15, 180, 30):
                        pi = Point((xi, yi))
                        if pi.intersects(overlap_region):
                            ax.scatter(xi, yi, color='green', alpha=0.8, zorder=1)
                            trace_x = xi+np.array(disp_box[grism])*ps['x']
                            trace_y = np.array([yi,yi])
                            trace_x, trace_y = threedhst.utils.xyrot(trace_x, trace_y, pair[j], x0=xi, y0=yi)
                            #ax.plot(trace_x, trace_y, color='green', alpha=0.4, linewidth=2, zorder=1)
                            trace_line = LineString([(trace_x[0], trace_y[0]), (trace_x[1], trace_y[1])])
                            trace_in_full = trace_line.buffer(1, resolution=8)
                            ax.add_patch(PolygonPatch(trace_in_full, fc='black', ec='None', alpha=0.15, zorder=1))
                            trace_in_full = trace_line.buffer(1, resolution=8).intersection(full_poly)
                            ax.add_patch(PolygonPatch(trace_in_full, fc='green', ec='None', alpha=0.5, zorder=1))
                            
    unicorn.plotting.savefig(fig, 'orient_demo.pdf')

    #### Testing
    if False:
        line = LineString([(-2,0.5), (0.5, 0.5)])
        full_box_x = np.array([0, 1, 1, 0, 0])
        full_box_y = np.array([0, 0, 1, 1, 0])
        box = Polygon(np.array([full_box_x, full_box_y]).T)
        line_within = line.buffer(0.02, resolution=8).intersection(box)
        
        plt.figure()
        ax = plt.gca()
        ax.add_patch(PolygonPatch(box, fc='black', ec='black', alpha=0.2, zorder=-2))
        ax.add_patch(PolygonPatch(line_within, fc='green', ec='None', alpha=0.2, zorder=-1))
        plt.xlim(-2,2)
        plt.ylim(-2,2)
        
        
def ISR_figures_coverage():
    """
    Show overlap coverage with/without shifting
    """
    
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    
    thetas = np.arange(-180, 181, 2)
    area_aper = thetas*0.
    area_cov = thetas*0.
    
    f_spec = 0.5

    fig = unicorn.plotting.plot_init(xs=8, aspect=1./2, square=True, left=0.09, wspace=0., bottom=0.1, right=0.01, top=0.01, use_tex=True, NO_GUI=True)
    
    ax, grism = fig.add_subplot(121), 'G102'
    ax, grism = fig.add_subplot(122), 'G141'
    
    for i in thetas:
        result = mywfc3.orient.overlap(theta=thetas[i], grism=grism, have_mosaic=True, plot=False, aper='GRISM1024', f_spec=f_spec)
        area_aper[i], area_cov[i] = result[2], result[5] 
        
    ax.plot(thetas, area_aper, label='GRISM1024', color='red')
    ax.plot(thetas, area_cov, label='Recenter', color='blue')
    
    for i in thetas:
        result = mywfc3.orient.overlap(theta=thetas[i], grism=grism, have_mosaic=False, plot=False, aper='GRISM1024', f_spec=f_spec)
        area_aper[i], area_cov[i] = result[2], result[5] 
        
    ax.plot(thetas, area_aper, label='GRISM1024, no edge' , color='red', linestyle='--')
    ax.plot(thetas, area_cov, label='Recenter, no edge', color='blue', linestyle='--')
    ax.legend(loc='lower center', title=grism + r', $f > %0.2f$' %(f_spec), fontsize=9)

    axes = fig.axes
    axes[0].set_ylabel('Overlap fraction')
    axes[1].set_yticklabels([])
    
    for ax in axes:
        ax.set_xlabel(r'$\theta = \Delta$ ORIENT (deg)')
        ax.xaxis.set_major_locator(MultipleLocator(45))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid()
        ax.set_ylim(0.55,1.01)
    
    unicorn.plotting.savefig(fig, 'orient_area_fraction.pdf'); plt.close()
    
    
