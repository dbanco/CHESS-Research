# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:50:59 2024

@author: dpqb1
"""
import numpy as np

def Ax_ft(A0, x):
    """
    GPU version: Computes forward convolution Ax = sum_k ifftn(A0[...,k] * fftn(x[...,k]))
    Supports both 1D (N x K) and 2D (N1 x N2 x K) data.

    Inputs:
    - A0: (..., K) np.array
    - x:  (..., K) np.array

    Output:
    - Ax: (...,) np.array
    """
    Ax = np.zeros(A0.shape[:-1], dtype=np.float64)
    x_ft = np.fft.fftn(x, axes=range(x.ndim - 1))

    for k in range(x.shape[-1]):
        Ax += np.real(np.fft.ifftn(A0[..., k] * x_ft[..., k], axes=range(x.ndim - 1)))

    return Ax

def AtR_ft(A0, R):
    """
    GPU version: Computes transpose operation AtR = [ifftn(conj(A0[...,k]) * fftn(R)) for k]
    Supports both 1D (N,) and 2D (N1 x N2) residuals.

    Inputs:
    - A0: (..., K) np.array
    - R:  (...,) np.array

    Output:
    - AtR: (..., K) np.array
    """
    AtR = np.zeros_like(A0, dtype=np.complex128)
    R_ft = np.fft.fftn(R, axes=range(R.ndim))

    for k in range(A0.shape[-1]):
        AtR[..., k] = np.fft.ifftn(np.conj(A0[..., k]) * R_ft, axes=range(R.ndim))

    return np.real(AtR)

def Ax_ft_Time(A0, X, Bnorms=None):
    """
    Compute A x where A is a stack of Fourier basis filters.
    Works for both 1D and 2D inputs.

    Inputs:
    - A0: N x K or N1 x N2 x K (fft of dictionary)
    - X: N x K x T or N1 x N2 x K x T
    - Bnorms: T-length array

    Output:
    - Y: N x T or N1 x N2 x T
    """
    T = X.shape[-1]
    is_2D = X.ndim == 4

    if is_2D:
        N1, N2, K = A0.shape
        Y = np.zeros((N1, N2, T))

        for t in range(T):
            norm = Bnorms[t] if Bnorms is not None else 1.0
            y_t = np.zeros((N1, N2), dtype=complex)

            for k in range(K):
                xk = np.fft.fft2(X[:, :, k, t])
                y_t += A0[:, :, k] * xk / norm

            Y[:, :, t] = np.real(np.fft.ifft2(y_t))

    else:
        N, K = A0.shape
        Y = np.zeros((N, T))

        for t in range(T):
            norm = Bnorms[t] if Bnorms is not None else 1.0
            y_t = np.zeros(N, dtype=complex)

            for k in range(K):
                xk = np.fft.fft(X[:, k, t])
                y_t += A0[:, k] * xk / norm

            Y[:, t] = np.real(np.fft.ifft(y_t))

    return Y


def AtB_ft_Time(A0, B, Bnorms):
    """
    Compute A^T B for both 1D and 2D.

    Inputs:
    - A0: N x K or N1 x N2 x K
    - B: N x T or N1 x N2 x T
    - Bnorms: T-length array

    Output:
    - AtB: N x K x T or N1 x N2 x K x T
    """
    T = B.shape[-1]
    is_2D = B.ndim == 3

    if is_2D:
        N1, N2, K = A0.shape
        AtB = np.zeros((N1, N2, K, T), dtype=complex)

        for t in range(T):
            norm = Bnorms[t] if Bnorms is not None else 1.0
            Bt_ft = np.fft.fft2(B[:, :, t]) / norm

            for k in range(K):
                AtB[:, :, k, t] = np.fft.ifft2(np.conj(A0[:, :, k]) * Bt_ft)

        return np.real(AtB)

    else:
        N, K = A0.shape
        AtB = np.zeros((N, K, T), dtype=complex)

        for t in range(T):
            norm = Bnorms[t] if Bnorms is not None else 1.0
            Bt_ft = np.fft.fft(B[:, t]) / norm

            for k in range(K):
                AtB[:, k, t] = np.fft.ifft(np.conj(A0[:, k]) * Bt_ft)

        return np.real(AtB)

def DiffPhiX(X, N=None, K=None, T=None):
    """
    Apply sum over spatial dimensions and temporal difference.
    Works for both 1D (N x K x T) and 2D (N1 x N2 x K x T).

    Output:
    - DiffPhi: K x (T-1)
    """

    dims = X.ndim
    if dims == 3:
        # 1D case: N x K x T
        N, K, T = X.shape if N is None else (N, K, T)
        PhiX = np.sum(X, axis=0)  # sum over N
    elif dims == 4:
        # 2D case: N1 x N2 x K x T
        N1, N2, K, T = X.shape
        PhiX = np.sum(X, axis=(0, 1))  # sum over N1, N2
    else:
        raise ValueError("Unsupported input dimensionality for DiffPhiX.")

    # Temporal difference
    DiffPhi = PhiX[:, 1:] - PhiX[:, :-1]  # shape: K x (T-1)

    return DiffPhi

def PhiTranDiffTran(R, shape):
    """
    Adjoint of temporal difference and summing operator Phi.
    Works for both 1D and 2D.

    Inputs:
    - R: K x (T-1) array (temporal residuals)
    - shape: Tuple (N, K, T) or (N1, N2, K, T)

    Output:
    - PtDtR: shape-matched tensor with adjoint operation applied
    """

    if len(shape) == 3:
        N, K, T = shape
        spatial_shape = (N,)
    elif len(shape) == 4:
        N1, N2, K, T = shape
        spatial_shape = (N1, N2)
    else:
        raise ValueError("Shape must be 3D or 4D for PhiTranDiffTran.")

    # Compute D^T R
    DtR = np.zeros((K, T))
    DtR[:, 0] = -R[:, 0]
    DtR[:, -1] = R[:, -1]
    DtR[:, 1:-1] = R[:, :-1] - R[:, 1:]

    # Expand back across spatial domain
    full_shape = spatial_shape + (K, T)
    PtDtR = np.tile(DtR[np.newaxis, ...], spatial_shape + (1,) * 2)

    return PtDtR