
from __future__ import division

import numpy as np
cimport numpy as np

DTYPE = np.double
ITYPE = np.int64

ctypedef np.double_t DTYPE_t
ctypedef np.uint_t UINT_t
ctypedef np.int_t INT_t
ctypedef np.int64_t LINT_t

import cython

cdef extern from "math.h":
    double sqrt(double x)
    double exp(double x)
    

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.embedsignature(True)
def disperse_grism_object(np.ndarray[DTYPE_t, ndim=2] flam, np.ndarray[LINT_t, ndim=1] idxl, np.ndarray[DTYPE_t, ndim=1] yfrac, np.ndarray[DTYPE_t, ndim=1] ysens, np.ndarray[DTYPE_t, ndim=1] full, np.ndarray[LINT_t, ndim=1] x0, np.ndarray[LINT_t, ndim=1] sh):

    cdef int i,j
    cdef unsigned int nk,k
    cdef double fl_ij
    
    nk = len(idxl)
        
    for i in range(0-sh[1], sh[1]):
        for j in range(0-sh[0], sh[0]):
            fl_ij = flam[x0[1]+j, x0[0]+i]/1.e-17
            for k in range(nk):
                full[idxl[k]+j*1014+i] += ysens[k]*fl_ij*yfrac[k]
                full[idxl[k]+(j-1)*1014+i] += ysens[k]*fl_ij*(1-yfrac[k])
    
    return True
    