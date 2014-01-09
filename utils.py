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
    
     