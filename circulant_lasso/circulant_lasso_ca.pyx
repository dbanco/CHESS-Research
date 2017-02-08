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


cdef inline floating fsign(floating f) nogil:
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
                 np.ndarray[np.complex128, ndim=2, mode='fortran'] Aft,
                 np.ndarray[floating, ndim=2, mode='fortran'] A_tran_ft,
                 np.ndarray[int, ndim=1] maxima,
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
       
#    assert A.dtype == floating and b.dtype == floating    
    
    cdef int m = A.shape[0]
    cdef int n = A.shape[1]
    cdef int num_vars = Aft.shape[1]
    cdef int num_maxima = maxima.shape[0]

    # precompute sum of sqaures of columns
    cdef np.ndarray[floating, ndim=1] z = np.sum(A**2,0)
   
    # get the number of tasks indirectly, using strides
    cdef int n_tasks = b.strides[0] / sizeof(floating)   
   
    # initialize x
    cdef np.ndarray[floating, ndim=1] x = np.zeros((n), dtype=dtype)
    cdef np.ndarray[floating, ndim=1] x_old = np.zeros((num_maxima), dtype=dtype)       
    cdef np.ndarray[floating, ndim=1] x_pad = np.zeros((m), dtype=dtype)        
    cdef np.ndarray[floating, ndim=1] AtR = np.zeros((n), dtype=dtype)
    cdef np.ndarray[floating, ndim=1] R = np.zeros((m), dtype=dtype)
    
    cdef floating *b_data = <floating*> b.data
    cdef floating *x_data = <floating*> x.data
    cdef floating *x_old_data = <floating*> x.data
    cdef floating *R_data = <floating*> R.data 
    cdef floating *AtR_data = <floating*> AtR.data
    
    cdef floating *A_data = <floating*> A.data
    cdef np.complex128 *Aft_data = <floating*> Aft.data 
    cdef np.complex128 *A_tran_ft_data = <floating*> A_tran_ft.data 
    
    cdef np.ndarray[np.complex128, ndim=1] x_ft = np.zeros((m), dtype=dtype)
    cdef np.ndarray[np.complex128, ndim=1] Rft = np.zeros((m), dtype=dtype)       
       
    cdef np.complex128 *x_ft_data = <floating*> x_ft.data
    cdef np.complex128 *Rft_data = <floating*> Rft.data    
    
    cdef floating rho    
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
       
    # Constant variables
    cdef floating neg_one = -1
    cdef floating zero_f = 0.0
    cdef floating one_f = 1.0
    cdef int zero = 0
    cdef int one = 1            
       
    # R = b - np.dot(A,x) but x is zeros   
    for i in range(m):
        R[i] = b[i]
        Rft[i] = b[i]
    
    # tol *= np.dot(b,b)
    tol *= dot(&m, b_data, &n_tasks, b_data, &n_tasks)
 
    for iters in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0
        
        for ii in range(num_vars):
            x_old = x[ii*num_maxima:(ii+1)*num_maxima]
            x_pad = np.zeros((m), dtype=dtype)
            x_pad[maxima] = x_old
            x_ft = np.fft.fft(x_pad)
            
            # Remove contribution to residual of current coefficient
            # R += x_ii*A[:,ii]
            sbmv("L", &m, &zero, &one_f, &Aft_data[ii * m], &one, x_ft_data, &one, &zero_f, Rft_data, &one);
            R = np.multiply            
                        
            R += np.fft.ifft(Rft)
            
            # Soft threshold update with optional positivity constraint
            for iii in range(num_maxima):   
                # rho = np.dot(A[:,ii],R)
                rho = dot(&m, &A_data[(iii+ii*num_maxima) * m], &one, R_data, &one )            
            
                if positive and rho < 0:
                    x[iii+num_maxima*ii] = 0.0
                else:
                    x[iii+num_maxima*ii] = (fsign(rho) * fmax(fabs(rho) - alpha, 0)/ (z[ii]))
            
            # Add contribution to residual of updated coefficient
            # R -= x[ii]*A[:,ii]
            x_pad = np.zeros((m), dtype=dtype)
            x_pad[maxima] = x[ii*num_maxima:(ii+1)*num_maxima]
            x_ft = np.fft.fft(x_pad)
            
            sbmv("L", &m, &zero, &one_f, &Aft_data[ii * m], &one, x_ft_data, &one, &zero_f, Rft_data, &one);
            R -= np.fft.ifft(Rft)
            
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
            for i in range(n):
                AtR[i] = dot(&m, &A_data[i * m], &one, R_data, &one)

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
    