"""
Utility scripts
"""

import glob
import os

def gzfile(filename):
    """
    Append '.gz' to filename if gzipped version of the file exists (and not gzipped doesn't)
    """
    if os.path.exists(filename):
        return filename
    
    for ext in ['.gz', '.Z']:
        if os.path.exists(filename+ext):
            return filename+ext
            
    return None
    
def wfc3_history():
    """
    Read a CSV file of all WFC3 observations searched from the 
    archive (as of Jan 2, 2014).
    """
    from astropy.io import ascii
    
    catalog = os.path.join(os.path.dirname(__file__), 'data/hst_search.txt.gz')
    print 'Read %s' %(catalog)
    tab = ascii.read(catalog, data_start=2, delimiter=',')
    tab.rename_column('RA (J2000)', 'ra')
    tab.rename_column('Dec (J2000)', 'dec')
    tab.rename_column('Filters/Gratings', 'filter')
    
    return tab
    
def parse_history():
    import numpy as np
    import mywfc3.utils
    t = mywfc3.utils.wfc3_history()
    
    progs = np.unique(t['Proposal ID'])
    programs = {}
    for prog in progs:
        mat = t['Proposal ID'] == prog
        filters = list(np.unique(t['filter'][mat]))
        f = {}
        for filter in filters:
            mfilt = mat & (t['filter'] == filter)
            f[filter] = np.round(np.sum(t['Exp Time'][mfilt])/5500.*2)
        #
        programs[prog] = {'PI':t['PI Last Name'][mat][0], 'NORB':np.round(np.sum(t['Exp Time'][mat])/5500.*2), 'filters':f}
        
    #### Grism programs
    for prog in progs:
        p = programs[prog]
        if 'G141' in p['filters'].keys():
            if p['filters']['G141'] > 5:
                print prog, p
    
    ## F105 / F110W / F098M
    for prog in progs:
        p = programs[prog]
        if 'F105W' in p['filters'].keys():
            if p['filters']['F105W'] > 2:
                print prog, p

def esa_fetch(program=14227):
    """
    Fetch data from the ESA HST archive
    """
    
    url = "http://archives.esac.esa.int/ehst-sl-server/servlet/metadata-action?RESOURCE_CLASS=OBSERVATION&PROPOSAL.PROPOSAL.PROPOSAL_ID=%d&SELECTED_FIELDS=OBSERVATION&RETURN_TYPE=CSV" %(program)
    os.system('curl -o esa.result.csv "%s"' %(url))
    
def fetch_nearby(radec=[40.055709,-1.6262032], r=3, instruments=['WFC3'], min_exptime=50):
    
    from threedhst import catIO
    
    dd = r/60.
    dr = dd/np.cos(radec[1]/180*np.pi)
    
    instrument_query = ' OR '.join(['INSTRUMENT.INSTRUMENT_NAME = \'%s\'' %(inst) for inst in instruments])
    
    url = r"http://archives.esac.esa.int/ehst-sl-server/servlet/metadata-action?RESOURCE_CLASS=SIMPLE_OBSERVATION&SELECTED_FIELDS=SIMPLE_OBSERVATION,INSTRUMENT.INSTRUMENT_NAME,OPTICAL_ELEMENT.OPTICAL_ELEMENT_NAME,PROPOSAL.PROPOSAL_ID,TARGET.TARGET_NAME,&QUERY=(SIMPLE_OBSERVATION.INTENT = 'Science' AND (%s) AND (POSITION.DEC > %.4f AND POSITION.DEC < %.4f AND POSITION.RA > %.4f AND POSITION.RA < %.4f) AND (SIMPLE_OBSERVATION.EXPOSURE_DURATION > %.1f))&ORDER_BY=SIMPLE_OBSERVATION.START_TIME&RETURN_TYPE=CSV" %(instrument_query, radec[1]-dd, radec[1]+dd, radec[0]-dr, radec[0]+dr, min_exptime)
        
    os.system('curl -o esa.csv "%s"' %(url.replace(' ', '%20')))
        
    query = catIO.Table('esa.csv')
    
    ### Get unique visits
    root = {'WFC3':'I', 'ACS':'J'}
    
    base_names = np.unique(['%s%s%s' %(root[l['INSTRUMENT_NAME']], l['PROGRAM_ID'], l['SET_ID']) for l in query])
    
    for base in base_names:
        for ext in [1,2,3,4,5,6]:
            for extn in ['ASN', 'JIT', 'JIF']:
                file_base = "%s0%d0_%s" %(base, ext, extn)
                file_url = r"http://archives.esac.esa.int/ehst-sl-server/servlet/data-action?ARTIFACT_ID=%s" %(file_base)
                os.system('curl -o %s.fits.gz "%s"' %(file_base.lower(), file_url))
            
def program_regions(program=13420, query=None):
    """
    Make region files from FLTs in the quicklook archive
    """
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    import glob
    import os
    
    if query is not None:
        files=glob.glob(query)
    else:
        files=glob.glob('/grp/hst/wfc3a/GO_Links/%d/Visit*/*flt.fits' %(program))
    
    files.sort()
    fp = open('%s.reg' %(program),'w')
    fp.write('fk5\n')
    
    old_visit = '--'
    for file in files:
        print file
        im = pyfits.open(file)

        visit = os.path.basename(file)[4:6]
        if visit == old_visit:
            visit = ''
        else:
            old_visit = visit

        if im[0].header['DETECTOR'] == 'IR':
            sci_ext = [1]
            if im[0].header['FILTER'] in ['G102', 'G141']:
                color='cyan'
            else:
                color='green'
        else:
            sci_ext = [1,4]
            color='magenta'
            
        for ext in sci_ext:
            wcs = pywcs.WCS(im[ext].header)
            footprint = wcs.calc_footprint(undistort=True)
            poly_coord = ['%.6f, %.6f' %(rd[0], rd[1]) for rd in footprint]
            fp.write('polygon(%s) # color=%s width=2 text={%s}\n' %(','.join(poly_coord), color, visit))
        
        ## Death-star for IR
        if len(sci_ext) == 1:
            rd = wcs.wcs_pix2world([357.8],[54.7],1)
            fp.write('circle(%.6f,%.6f,2.9") # color=%s\n' %(rd[0][0], rd[1][0], color))
            
    fp.close()
     