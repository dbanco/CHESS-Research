# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 16:17:49 2016

@author: Dan
"""


from libc.math cimport fabs
cimport numpy as np
import numpy as np
import numpy.linalg as linalg

cimport cython
from cpython cimport bool
from cython cimport floating
import warnings

ctypedef np.float64_t DOUBLE
ctypedef np.uint32_t UINT32_t
ctypedef floating (*DOT)(int N, floating *X, int incX, floating *Y,
                         int incY) nogil
ctypedef void (*AXPY)(int N, floating alpha, floating *X, int incX,
                      floating *Y, int incY) nogil
ctypedef floating (*ASUM)(int N, floating *X, int incX) nogil

np.import_array()

# The following two functions are shamelessly copied from the tree code.

cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    RAND_R_MAX = 0x7FFFFFFF


cdef inline UINT32_t our_rand_r(UINT32_t* seed) nogil:
    seed[0] ^= <UINT32_t>(seed[0] << 13)
    seed[0] ^= <UINT32_t>(seed[0] >> 17)
    seed[0] ^= <UINT32_t>(seed[0] << 5)

    return seed[0] % (<UINT32_t>RAND_R_MAX + 1)


cdef inline UINT32_t rand_int(UINT32_t end, UINT32_t* random_state) nogil:
    """Generate a random integer in [0; end)."""
    return our_rand_r(random_state) % end


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
def lasso_coordinate_descent(floating alpha,
                             np.ndarray[floating, ndim=2, mode='fortran'] X,
                             np.ndarray[floating, ndim=1, mode='c'] y,
                             int max_iter, floating tol,
                             object rng, bint random=0, bint positive=0):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression
        We minimize
        (1/2) * norm(y - X w, 2)^2 + alpha norm(w, 1)
    """

    # fused types version of BLAS functions
    cdef DOT dot
    cdef AXPY axpy
    cdef ASUM asum

    if floating is float:
        dtype = np.float32
        dot = sdot
        axpy = saxpy
        asum = sasum
    else:
        dtype = np.float64
        dot = ddot
        axpy = daxpy
        asum = dasum

    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    # get the number of tasks indirectly, using strides
    cdef unsigned int n_tasks = y.strides[0] / sizeof(floating)

    # compute norms of the columns of X
    cdef np.ndarray[floating, ndim=1] norm_cols_X = (X**2).sum(axis=0)
    cdef np.ndarray[floating, ndim=1] w = np.empty(n_features, dtype=dtype)
    # initial value of the residuals
    cdef np.ndarray[floating, ndim=1] R = np.empty(n_samples, dtype=dtype)
    cdef np.ndarray[floating, ndim=1] XtA = np.empty(n_features, dtype=dtype)

    cdef floating tmp
    cdef floating w_ii
    cdef floating d_w_max
    cdef floating w_max
    cdef floating d_w_ii
    cdef floating gap = tol + 1.0
    cdef floating d_w_tol = tol
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating w_norm2
    cdef floating l1_norm
    cdef floating const
    cdef floating A_norm2
    cdef unsigned int ii
    cdef unsigned int i
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed

    cdef floating *X_data = <floating*> X.data
    cdef floating *y_data = <floating*> y.data
    cdef floating *w_data = <floating*> w.data
    cdef floating *R_data = <floating*> R.data
    cdef floating *XtA_data = <floating*> XtA.data

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

    with nogil:
        # R = y - np.dot(X, w)
        for i in range(n_samples):
            R[i] = y[i] - dot(n_features, &X_data[i], n_samples, w_data, 1)

        # tol *= np.dot(y, y)
        tol *= dot(n_samples, y_data, n_tasks, y_data, n_tasks)

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm_cols_X[ii] == 0.0:
                    continue

                w_ii = w[ii]  # Store previous value

                if w_ii != 0.0:
                    # R += w_ii * X[:,ii]
                    axpy(n_samples, w_ii, &X_data[ii * n_samples], 1,
                         R_data, 1)

                # tmp = (X[:,ii]*R).sum()
                tmp = dot(n_samples, &X_data[ii * n_samples], 1, R_data, 1)

                if positive and tmp < 0:
                    w[ii] = 0.0
                else:
                    w[ii] = (fsign(tmp) * fmax(fabs(tmp) - alpha, 0)
                             / (norm_cols_X[ii]))

                if w[ii] != 0.0:
                    # R -=  w[ii] * X[:,ii] # Update residual
                    axpy(n_samples, -w[ii], &X_data[ii * n_samples], 1,
                         R_data, 1)

                # update the maximum absolute coefficient update
                d_w_ii = fabs(w[ii] - w_ii)
                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                if fabs(w[ii]) > w_max:
                    w_max = fabs(w[ii])

            if (w_max == 0.0 or
                d_w_max / w_max < d_w_tol or
                n_iter == max_iter - 1):
                # the biggest coordinate update of this iteration was smaller
                # than the tolerance: check the duality gap as ultimate
                # stopping criterion

                # XtA = np.dot(X.T, R)
                for i in range(n_features):
                    XtA[i] = dot(n_samples, &X_data[i * n_samples],
                                 1, R_data, 1)

                if positive:
                    dual_norm_XtA = max(n_features, XtA_data)
                else:
                    dual_norm_XtA = abs_max(n_features, XtA_data)

                # R_norm2 = np.dot(R, R)
                R_norm2 = dot(n_samples, R_data, 1, R_data, 1)

                # w_norm2 = np.dot(w, w)
                w_norm2 = dot(n_features, w_data, 1, w_data, 1)

                if (dual_norm_XtA > alpha):
                    const = alpha / dual_norm_XtA
                    A_norm2 = R_norm2 * (const ** 2)
                    gap = 0.5 * (R_norm2 + A_norm2)
                else:
                    const = 1.0
                    gap = R_norm2

                l1_norm = asum(n_features, w_data, 1)

                # np.dot(R.T, y)
                gap += (alpha * l1_norm
                        - const * dot(n_samples, R_data, 1, y_data, n_tasks))

                if gap < tol:
                    # return if we reached desired tolerance
                    break
    return w, gap, tol, n_iter + 1