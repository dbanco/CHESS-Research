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

######### NEED TO REIMPLEMENT? ######### 
"Soft thresholding function with vector input and output"
def soft_thresh(x, l):
    return np.sign(x) * np.maximum(np.abs(x) - l, 0.) 
    
"Vector formatting helper function for circulant matrix-vector product"
def x_to_x_ft(x, maxima, m, n):
    """ 
    Inputs:
        x       coefficient vector
        maxima  vector of locations of maxima in signal
        m       number of rows in x_ft
        n       number of columns in x_ft
    Output:
        x_ft    Each column is the fourier transform of the coefficients
                corresponding to the circulant matrix Ai
    """
    x_pad = np.zeros((m,n))
    x_ft = np.zeros((m,n))
    num_maxima = len(maxima)
    for i in range(n):
        x_pad[maxima,i] = x[i*num_maxima:(i+1)*num_maxima]      
        
    x_ft = np.fft.fft(x_pad,axis=0)
    
    return x_ft

"Circulant matrix-vector product subroutine"
def Ax_ft(A0ft, x, maxima, m, n):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        x       Coefficient vector
        maxima  vector of locations of maxima in signal
        m       number of rows in x_ft
        n       number of columns in x_ft
    """
    Ax = np.zeros((A0ft.shape[0]))
    
    x_ft = x_to_x_ft(x,maxima,m,n)
    
    for ii in range(A0ft.shape[1]):
        Ax += np.fft.ifft(np.multiply(A0ft[:,ii],x_ft[:,ii])).real
        
    return Ax

"Transposed circulant matrix-vector product subroutine"
def ATb_ft(A0_tran_ft,R,num_var,maxima):
    # Inputs:
    #    A0_tran_ft    Each column is the first row of a circulant matrix Ai
    #    R             Residual vector 
    # Output:
    #    ATb           Remember that this will be zero padded

    num_maxima = len(maxima)
    ATb = np.zeros((num_var*num_maxima))
    R_ft = np.fft.fft(R)          
    
    # num_vars    
    for ii in range(num_var):
        ATb[ii*num_maxima:(ii+1)*num_maxima] = np.fft.ifft(np.multiply(A0_tran_ft[:,ii],R_ft)).real[maxima]   

    return ATb
#########                      #########

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def circulant_fista(np.ndarray[floating, ndim=2, mode='fortran'] A, 
                    np.ndarray[floating, ndim=1] b,
                     floating alpha, 
                     int max_iter,
                     bint positive=0):
    
    cdef DOT dot
    cdef AXPY axpy
    cdef ASUM asum
    if floating is float:
        dtype = np.float32
        dot = cy_blas.sdot
        axpy = cy_blas.saxpy
        asum = cy_blas.sasum
    else:
        dtype = np.float64
        dot = cy_blas.ddot
        axpy = cy_blas.daxpy
        asum = cy_blas.dasum
       
#    assert A.dtype == floating and b.dtype == floating    
    
    cdef int m = A.shape[0]
    cdef int n = A.shape[1]
    cdef int o = 1

    # precompute sum of sqaures of columns
    cdef np.ndarray[floating, ndim=1] z = np.sum(A**2,0)
   
    # get the number of tasks indirectly, using strides
    cdef int n_tasks = b.strides[0] / sizeof(floating)   
   
    # initialize x
    cdef np.ndarray[floating, ndim=1] x = np.zeros((n), dtype=dtype)
    cdef np.ndarray[floating, ndim=1] R = np.zeros((m), dtype=dtype)
    cdef np.ndarray[floating, ndim=1] AtR = np.zeros((n), dtype=dtype)
    
    cdef floating *A_data = <floating*> A.data
    cdef floating *b_data = <floating*> b.data
    cdef floating *x_data = <floating*> x.data
    cdef floating *R_data = <floating*> R.data 
    cdef floating *AtR_data = <floating*> AtR.data
    
    cdef floating rho    
    cdef floating x_max
    cdef floating d_x_max  
    cdef floating x_ii
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
    cdef floating neg_tmp
    # loop variables    
    cdef unsigned int iters = 0
    cdef unsigned int i = 0    
    cdef unsigned int ii = 0
       
    # R = b - np.dot(A,x)
    for i in range(m):
        R[i] = b[i] - dot(&n, &A_data[i], &m, x_data, &o)
    
    # tol *= np.dot(b,b)
    tol *= dot(&m, b_data, &n_tasks, b_data, &n_tasks)
 
    for iters in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0

        if z[ii] == 0.0:
            continue
        
        for ii in range(n):
            x_ii = x[ii]
            
            # Remove contribution to residual of current coefficient
            # R += x_ii*A[:,ii]
            if(x_ii != 0.0):
                axpy(&m, &x_ii, &A_data[ii * m], &o, R_data, &o)

            # rho = np.dot(A[:,ii],R)
            rho = dot(&m, &A_data[ii * m], &o, R_data, &o )
            
            # Soft threshold update with optional positivity constraint
            if positive and rho < 0:
                x[ii] = 0.0
            else:
                x[ii] = (fsign(rho) * fmax(fabs(rho) - alpha, 0)/ (z[ii]))
            
            # Add contribution to residual of updated coefficient
            # R -= x[ii]*A[:,ii]
            if(x[ii] != 0.0):
                neg_tmp = -x[ii]
                axpy(&m, &neg_tmp, &A_data[ii * m], &o, R_data, &o)
                
            # update the maximum absolute coefficient update
            d_x_ii = fabs(x[ii] - x_ii)
            if d_x_ii > d_x_max:
                d_x_max = d_x_ii

            if fabs(x[ii]) > x_max:
                x_max = fabs(x[ii])

            if (x_max == 0.0 or
                d_x_max / x_max < d_x_tol or
                iters == max_iter - 1):
                # the biggest coordinate update of this iteration was smaller
                # than the tolerance: check the duality gap as ultimate
                # stopping criterion

                # AtR = np.dot(A.T, R)
                for i in range(n):
                    AtR[i] = dot(&m, &A_data[i * m], &o, R_data, &o)

                if positive:
                    dual_norm_AtR = max(n, AtR_data)
                else:
                    dual_norm_AtR = abs_max(n, AtR_data)

                # R_norm2 = np.dot(R, R)
                R_norm2 = dot(&m, R_data, &o, R_data, &o)

                # w_norm2 = np.dot(w, w)
                x_norm2 = dot(&n, x_data, &o, x_data, &o)

                if (dual_norm_AtR > alpha):
                    const = alpha / dual_norm_AtR
                    A_norm2 = R_norm2 * (const ** 2)
                    gap = 0.5 * (R_norm2 + A_norm2)
                else:
                    const = 1.0
                    gap = R_norm2

                l1_norm = asum(&n, x_data, &o)

                # np.dot(R.T, y)
                gap += (alpha * l1_norm
                        - const * dot(&m, R_data, &o, b_data, &n_tasks))

                if gap < tol:
                    # return if we reached desired tolerance
                    break
    return x
    
    cdef floating start_time
    if benchmark: 
        start_time = time.time()
            
    cdef int m = A.shape[0]    
    cdef int n = A0.shape[1]
    
    cdef np.ndarray[floating, ndim=2, mode='fortran'] A0ft = np.fft.fft(A0,axis=0)    
    cdef np.ndarray[floating, ndim=2, mode='fortran'] A0_tran_ft = np.fft.fft(A0_tran,axis=0)  
    
    x = np.zeros(A.shape[1])
    t = 1
    z = x.copy()
    L = norm(A)**2
    for it in range(maxit):
        xold = x.copy()
        
        # Arrange x coefficents as matrix in fourier domain 
        R = b - Ax_ft(A0ft, z, maxima, m, n)
        z = z + ATb_ft(A0_tran_ft,R,n,maxima)/L
        if positive:
                z[z<0]=0
        #z = z + A.T.dot(R)/L
        x = soft_thresh(z,1/L)
        t0 = t
        t = (1. + sqrt(1. + 4. * t ** 2)) / 2.
        z = x + ((t0 - 1.) / t) * (x - xold)

    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x    
    
    