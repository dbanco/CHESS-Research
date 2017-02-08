from libc.math cimport fabs
cimport numpy as np
import numpy as np

from numpy.linalg import norm 
cimport scipy.linalg.cython_blas as cy_blas
import warnings

cimport cython
from cython cimport floating
from cpython cimport bool
ctypedef np.float64_t DOUBLE
ctypedef np.uint32_t UINT32_t
ctypedef floating (*DOT)(int *N, floating *X, int *incX, floating *Y,
                         int *incY) nogil
ctypedef void (*AXPY)(int *N, floating *alpha, floating *X, int *incX,
                      floating *Y, int *incY) nogil
ctypedef floating (*ASUM)(int *N, floating *X, int *incX) nogil
ctypedef void (*SBMV)(char *uplo, int *n, int *k, floating *alpha, floating *a, int *lda, floating *x, int *incx, floating *beta, floating *y, int *incy) nogil

np.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def element_multiply(np.ndarray[floating, ndim=1] a_fft_re, 
                     np.ndarray[floating, ndim=1] a_fft_im, 
                     np.ndarray[floating, ndim=1] b_fft_re,
                     np.ndarray[floating, ndim=1] b_fft_im):
    
    cdef DOT dot
    cdef AXPY axpy
    cdef ASUM asum
    cdef SBMV sbmv
    
    if floating is float:
        dtype = np.float32
        dot = cy_blas.sdot
        axpy = cy_blas.saxpy
        asum = cy_blas.sasum
        sbmv = cy_blas.ssbmv
    else:
        dtype = np.float64
        dot = cy_blas.ddot
        axpy = cy_blas.daxpy
        asum = cy_blas.dasum
        sbmv = cy_blas.dsbmv

    cdef int n = a_fft_re.shape[0]

    cdef floating *a_fft_re_data = <floating*> a_fft_re.data
    cdef floating *a_fft_im_data = <floating*> a_fft_im.data
    cdef floating *b_fft_re_data = <floating*> b_fft_re.data
    cdef floating *b_fft_im_data = <floating*> b_fft_im.data

    # Constant variables
    cdef floating neg_one = -1
    cdef floating zero_f = 0.0
    cdef floating one_f = 1.0
    cdef floating neg_one_f = -1.0
    cdef int zero = 0
    cdef int one = 1            
       
    # Elementwise multiply to get real part   
    cdef np.ndarray[floating, ndim=1] c_fft_re = np.zeros(n, dtype=dtype)
    cdef floating *c_fft_re_data = <floating*> c_fft_re.data    
    #sbmv("L", &n, &zero, &one_f,     a_fft_re_data, &one, b_fft_re_data, &one, &zero_f, c_fft_re_data, &one);
    #sbmv("L", &n, &zero, &neg_one_f, a_fft_im_data, &one, b_fft_im_data, &one, &one_f,  c_fft_re_data, &one);
    
    # Elementwise multiply to get imaginary part  
    cdef np.ndarray[floating, ndim=1] c_fft_im = np.zeros(n, dtype=dtype)
    cdef floating *c_fft_im_data = <floating*> c_fft_im.data
    #sbmv("L", &n, &zero, &one_f, a_fft_re_data, &one, b_fft_im_data, &one, &zero_f, c_fft_im_data, &one);
    #sbmv("L", &n, &zero, &one_f, a_fft_im_data, &one, b_fft_re_data, &one, &one_f,  c_fft_im_data, &one);
    
    for i in range(n):
        c_fft_re[i] += a_fft_re[i]*b_fft_re[i] - a_fft_im[i]*b_fft_im[i]
        c_fft_im[i] += a_fft_re[i]*b_fft_im[i] + a_fft_im[i]*b_fft_re[i]
    
    return c_fft_re, c_fft_im
    