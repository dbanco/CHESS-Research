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

cdef inline floating fmax(floating x, floating y) nogil:
    if x > y:
        return x
    return y


cdef floating fsign(floating f) nogil:
    if f == 0:
        return 0
    elif f > 0:
        return 1.0
    else:
        return -1.0


cdef floating abs_max(int n, floating* a) nogil:
    """np.max(np.abs(a))"""
    cdef int i
    cdef floating m = fabs(a[0])
    cdef floating d
    for i in range(1, n):
        d = fabs(a[i])
        if d > m:
            m = d
    return m

cdef floating max(int n, floating* a) nogil:
    """np.max(a)"""
    cdef int i
    cdef floating m = a[0]
    cdef floating d
    for i in range(1, n):
        d = a[i]
        if d > m:
            m = d
    return m

cdef floating diff_abs_max(int n, floating* a, floating* b) nogil:
    """np.max(np.abs(a - b))"""
    cdef int i
    cdef floating m = fabs(a[0] - b[0])
    cdef floating d
    for i in range(1, n):
        d = fabs(a[i] - b[i])
        if d > m:
            m = d
    return m

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def lasso_solver(np.ndarray[floating, ndim=2, mode='fortran'] A,
                 np.ndarray[floating, ndim=2, mode='fortran'] A0,
                 np.ndarray[floating, ndim=2, mode='fortran'] A0_tran,
                 maxima,
                 np.ndarray[floating, ndim=1] b,
                 floating alpha, 
                 int max_iter,
                 bint positive=0):
    
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
   
    cdef int m = A.shape[0]
    cdef int n = A.shape[1]
    cdef int num_vars = A0.shape[1]
    cdef int num_maxima = maxima.shape[0]

    # precompute sum of sqaures of columns
    cdef np.ndarray[floating, ndim=1] z = np.sum(A**2,0)
   
    # get the number of tasks indirectly, using strides
    cdef int n_tasks = b.strides[0] / sizeof(floating)   
   
    # initialize x
    cdef np.ndarray[floating, ndim=1] x = np.zeros(n, dtype=dtype)
    cdef np.ndarray[floating, ndim=1] x_old = np.zeros(num_maxima, dtype=dtype)       
    cdef np.ndarray[floating, ndim=1] x_pad = np.zeros(m, dtype=dtype)        
    cdef np.ndarray[floating, ndim=1] AtR = np.zeros(m, dtype=dtype)
    cdef np.ndarray[floating, ndim=1] R = np.zeros(m, dtype=dtype) 
    Rft = np.zeros(m, dtype=np.complex128) 
    cdef floating *b_data = <floating*> b.data
    cdef floating *x_data = <floating*> x.data
    cdef floating *x_old_data = <floating*> x.data
    cdef floating *R_data = <floating*> R.data 
    cdef floating *AtR_data = <floating*> AtR.data
    
    cdef floating *A_data = <floating*> A.data
    
    cdef np.ndarray[floating, ndim=1] rho = np.zeros(m, dtype=dtype)   
    cdef floating rho_i
    cdef floating x_max
    cdef floating d_x_max  
    cdef floating tol = 1e-4
    cdef floating d_x_ii
    cdef floating gap = tol + 1.0
    cdef floating d_x_tol = tol
    cdef floating dual_norm_AtR
    cdef floating R_norm2
    cdef floating x_norm2
    cdef floating l1_norm
    cdef floating const
    cdef floating A_norm2    
    
    # loop variables    
    cdef unsigned int iters = 0
    cdef unsigned int i = 0    
    cdef unsigned int ii = 0
    cdef unsigned int iii = 0       
       
    # Constant variables
    cdef floating neg_one = -1
    cdef floating zero_f = 0.0
    cdef floating one_f = 1.0
    cdef int zero = 0
    cdef int one = 1            
       
    Aft = np.fft.fft(A0)
    A_tran_ft = np.fft.fft(A0_tran)       
       
    # R = b - np.dot(A,x) but x is zeros   
    for i in range(m):
        R[i] = b[i]
        Rft[i] = b[i]
    
    # tol *= np.dot(b,b)
    tol *= dot(&m, b_data, &n_tasks, b_data, &n_tasks)
 
    for iters in range(max_iter):
        x_max = 0.0
        d_x_max = 0.0
        
        for ii in range(num_vars):
            # 1)            
            x_old = x[ii*num_maxima:(ii+1)*num_maxima]
            x_pad = np.zeros((m), dtype=dtype)
            x_pad[maxima] = x_old
            x_ft = np.fft.fft(x_pad)
            
            # Remove contribution to residual of current coefficient
            # R += x_ii*A[:,ii]
            R += np.fft.ifft(np.multiply(Aft[:,ii],x_ft)).real
  
            # 2)
            # Soft threshold update with optional positivity constraint
            rho = np.fft.ifft(np.multiply(A_tran_ft[:,ii],np.fft.fft(R))).real          
            for iii,max_i in enumerate(maxima):   
                # rho = np.dot(A[:,ii],R)
                rho_i = rho[max_i]
                if positive and rho_i < 0:
                    x[iii+num_maxima*ii] = 0.0
                else:
                    x[iii+num_maxima*ii] = (fsign(rho_i) * fmax(fabs(rho_i) - alpha, 0)/ (z[ii]))
            
            # 3)
            x_pad = np.zeros((m), dtype=dtype)
            x_pad[maxima] = x[ii*num_maxima:(ii+1)*num_maxima]
            x_ft = np.fft.fft(x_pad)
            
            # Add contribution to residual of updated coefficient
            # R -= x[ii]*A[:,ii]
            R -= np.fft.ifft(np.multiply(Aft[:,ii],x_ft)).real
            
            # update the maximum absolute coefficient update
            d_x_ii = diff_abs_max(num_maxima, &x_data[ii*num_maxima], x_old_data)
            if d_x_ii > d_x_max:
                d_x_max = d_x_ii

            if abs_max(num_maxima, &x[ii*num_maxima]) > x_max:
                x_max = abs_max(num_maxima, &x[ii*num_maxima])

        if (x_max == 0.0 or
            d_x_max / x_max < d_x_tol or
            iters == max_iter - 1):
            # the biggest coordinate update of this iteration was smaller
            # than the tolerance: check the duality gap as ultimate
            # stopping criterion

            # AtR = np.dot(A.T, R)
            AtR = np.zeros((m))
            Rft = np.fft.fft(R)

            for i in range(num_vars):
                AtR += np.real(np.fft.ifft(np.multiply(A_tran_ft[:,i],Rft)))


            if positive:
                dual_norm_AtR = max(n, AtR_data)
            else:
                dual_norm_AtR = abs_max(n, AtR_data)

            # R_norm2 = np.dot(R, R)
            R_norm2 = dot(&m, R_data, &one, R_data, &one)

            # x_norm2 = np.dot(x, x)
            x_norm2 = dot(&n, x_data, &one, x_data, &one)

            if (dual_norm_AtR > alpha):
                const = alpha / dual_norm_AtR
                A_norm2 = R_norm2 * (const ** 2)
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            l1_norm = asum(&n, x_data, &one)

            # np.dot(R.T, y)
            gap += (alpha * l1_norm
                    - const * dot(&m, R_data, &one, b_data, &n_tasks))

            if gap < tol:
                # return if we reached desired tolerance
                break
    return x
    