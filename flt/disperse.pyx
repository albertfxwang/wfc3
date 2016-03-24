
from __future__ import division

import numpy as np
cimport numpy as np

DTYPE = np.double
ITYPE = np.int64

ctypedef np.double_t DTYPE_t
ctypedef np.uint_t UINT_t
ctypedef np.int_t INT_t
ctypedef np.int64_t LINT_t
ctypedef np.float32_t SEGM_t

import cython

cdef extern from "math.h":
    double sqrt(double x)
    double exp(double x)
    

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.embedsignature(True)
def disperse_grism_object(np.ndarray[DTYPE_t, ndim=2] flam, np.ndarray[SEGM_t, ndim=2] segm, int seg_id, np.ndarray[LINT_t, ndim=1] idxl, np.ndarray[DTYPE_t, ndim=1] yfrac, np.ndarray[DTYPE_t, ndim=1] ysens, np.ndarray[DTYPE_t, ndim=1] full, np.ndarray[LINT_t, ndim=1] x0, np.ndarray[LINT_t, ndim=1] shd, np.ndarray[LINT_t, ndim=1] shg):

    cdef int i,j,k1,k2
    cdef unsigned int nk,nl,k
    cdef double fl_ij
    
    nk = len(idxl)
    nl = len(full)
    
    for i in range(0-shd[1], shd[1]):
        for j in range(0-shd[0], shd[0]):
            fl_ij = flam[x0[0]+j, x0[1]+i]/1.e-17
            if (fl_ij == 0) | (segm[x0[0]+j, x0[1]+i] != seg_id):
                continue
                
            for k in range(nk):
                k1 = idxl[k]+j*shg[1]+i
                if (k1 >= 0) & (k1 < nl):
                    full[k1] += ysens[k]*fl_ij*yfrac[k]
                    
                k2 = idxl[k]+(j-1)*shg[1]+i
                if (k2 >= 0) & (k2 < nl):
                    full[k2] += ysens[k]*fl_ij*(1-yfrac[k])
    
    return True

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.embedsignature(True)
def seg_flux(np.ndarray[DTYPE_t, ndim=2] flam, np.ndarray[LINT_t, ndim=1] idxl, np.ndarray[DTYPE_t, ndim=1] yfrac, np.ndarray[DTYPE_t, ndim=1] ysens, np.ndarray[DTYPE_t, ndim=1] full, np.ndarray[LINT_t, ndim=1] x0, np.ndarray[LINT_t, ndim=1] shd, np.ndarray[LINT_t, ndim=1] shg):
    pass