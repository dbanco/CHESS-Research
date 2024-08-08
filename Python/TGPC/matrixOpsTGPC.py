# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:50:59 2024

@author: dpqb1
"""
import numpy as np

def Ax_ft_1D(A0, x):
    """
    Computes matrix-vector product between A and x.
    Element-wise multiplies each basis function of A0 with fft(x).

    Inputs:
    - A0: N x K array
    - x: N x K array

    Outputs:
    - Ax: N x 1 array
    """

    Ax = np.zeros(A0.shape[0])

    x_ft = np.fft.fft(x, axis=0)

    for k in range(x.shape[1]):
        Ax += np.real(np.fft.ifft(A0[:, k] * x_ft[:, k]))

    return Ax

def AtR_ft_1D(A0, R):
    """
    Computes matrix-vector product between A transposed and R.
    Element-wise multiplies each basis function of A0 with fft(R).

    Inputs:
    - A0: N x K array
    - R: N x 1 array

    Outputs:
    - AtR: N x K array
    """

    AtR = np.zeros_like(A0)

    R_ft = np.fft.fft(R)

    for k in range(A0.shape[1]):
        y = np.fft.ifft(np.conj(A0[:, k]) * R_ft)
        AtR[:, k] = y

    return AtR

def Ax_ft_1D_Time(A0, X, Bnorms=None):
    """
    Compute the product Ax where A is the convolution operator represented by A0
    with each frame of X along the third dimension, and normalize each frame
    by the corresponding value in Bnorms.

    Inputs:
    - A0: N x K FFT of dictionary
    - X: N x K x T tensor
    - Bnorms: T x 1 array of norm values for each frame in X

    Output:
    - Y: N x T tensor
    """
    T = X.shape[2]
    N, _ = A0.shape

    Y = np.zeros((N, T))

    if Bnorms is not None:
        for t in range(T):
            Y[:, t] = Ax_ft_1D(A0 / Bnorms[t], X[:, :, t])
    else:
        for t in range(T):
            Y[:, t] = Ax_ft_1D(A0, X[:, :, t])

    return Y

def AtB_ft_1D_Time(A0, B, Bnorms):
    """
    Compute the product AtB where A is the convolution operator represented by A0
    with each frame of B along the second dimension, and normalize each frame
    by the corresponding value in Bnorms.

    Inputs:
    - A0: N x K FFT of dictionary
    - B: N x T array
    - Bnorms: T x 1 array of norm values for each frame in B

    Output:
    - AtB: N x K x T tensor
    """
    N, K = A0.shape
    T = B.shape[1]
    AtB = np.zeros((N, K, T))

    for t in range(T):
        AtB[:, :, t] = AtR_ft_1D(A0 / Bnorms[t], B[:, t])

    return AtB

def DiffPhiX_1D(X, N=None, K=None, T=None):
    """
    Apply summing over space and temporal difference matrix.

    Inputs:
    - X: N x K x T array
    - N: Integer
    - K: Integer (optional)
    - T: Integer (optional)

    Outputs:
    - DiffPhi: N x (T-1) array
    """

    # Reshape X if it is a vector
    if N is None:
        N, K, T = X.shape
    else:
        X = X.reshape((N, K, T))

    PhiX = np.sum(X, axis=0)
    DiffPhi = PhiX[:, 1:] - PhiX[:, :-1]

    return DiffPhi

def PhiTranDiffTran_1D(R, N, K=None, T=None):
    """
    Compute residual for conjugate gradient that includes difference matrix.

    Inputs:
    - R: K x T array
    - N: Integer
    - K: Integer (optional)
    - T: Integer (optional)

    Outputs:
    - PtDtR: N x K x T array
    """

    # Reshape R if it is a vector
    if K is None:
        K, T = R.shape
        T += 1
        reshapeFlag = 0
    else:
        R = R.reshape((K, T))
        T += 1
        reshapeFlag = 1

    DtR = np.zeros((K, T))
    DtR[:, 0] = -R[:, 0]
    DtR[:, -1] = R[:, -1]
    DtR[:, 1:-2] = R[:, :-2] - R[:, 1:-1]

    DtR = DtR.reshape((1, K, T))
    PtDtR = np.tile(DtR, (N, 1, 1))

    if reshapeFlag:
        PtDtR = PtDtR.reshape((1,) + PtDtR.shape)

    return PtDtR