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
    