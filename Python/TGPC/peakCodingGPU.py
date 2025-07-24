# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:48:32 2024

@author: dpqb1
"""

import numpy as np
import cupy as cp
import matrixOpsTGPCGPU as mat

def conjGrad_TVphi_1D(A0, B, Bnorms, X_init, YV, ZU, params):
    """
    Solves linear system:
    At A x + rho1 I + rho2(Phit Phi X) = At b + rho1(Y-V) + rho2(Z-U)

    Inputs:
    - A0: N x K fft of dictionary
    - B: N x T data
    - Bnorms: T x 1 norm values of data
    - X_init: N x K x T initial solution
    - YV: N x K x T Y-V
    - ZU: K x T-1 Z-U
    - params: Dictionary containing parameters:
        - rho1: Dual variable 1
        - rho2: Dual variable 2
        - conjGradIter: Max number of conjugate gradient iterations
        - cgEpsilon: Stopping threshold

    Outputs:
    - Xk: N x K x T solution
    - cgIters: Final number of iterations
    """

    # ADMM penalty parameter
    rho1 = params['rho1']
    rho2 = params['rho2']
    N = A0.shape[0]

    # Coefficeint Vectors
    Xk = X_init.copy()

    # Target Vectors
    AtB = mat.AtB_ft_Time(A0, B, Bnorms)
    PtDtZ = mat.PhiTranDiffTran(ZU, Xk.shape)

    # Initial Residual
    Rk = AtB - AtAx(A0, Xk, Bnorms) + rho2 * PtDtZ - rho2 * PtDtDPx(Xk) + rho1 * YV - rho1 * Xk
    Pk = Rk.copy()

    for i in range(params['conjGradIter']):
        Apk = AtAx(A0, Pk, Bnorms) + rho2 * PtDtDPx(Pk) + rho1 * Pk
        RkRk = np.sum(Rk * Rk)
        alphak = RkRk / np.sum(Pk * Apk)
        Xk = Xk + alphak * Pk
        Rkp1 = Rk - alphak * Apk
        if np.linalg.norm(Rkp1) < params['cgEpsilon']:
            break
        betak = np.sum(Rkp1 * Rkp1) / RkRk
        Pk = Rkp1 + betak * Pk
        Rk = Rkp1

    cgIters = i + 1
    return Xk, cgIters

def conjGrad_TVphi_2D(A0, B, Bnorms, X_init, YV, ZU, params):
    """
    Solves linear system:
    At A x + rho1 I + rho2(Phit Phi X) = At b + rho1(Y-V) + rho2(Z-U)

    Inputs:
    - A0: N1 x N2 x K fft of dictionary
    - B: N1 x N2 x T data
    - Bnorms: T x 1 norm values of data
    - X_init: N1 x N2 x K x T initial solution
    - YV: N1 x N2 x K x T Y-V
    - ZU: K x T-1 Z-U
    - params: Dictionary containing parameters:
        - rho1: Dual variable 1
        - rho2: Dual variable 2
        - conjGradIter: Max number of conjugate gradient iterations
        - cgEpsilon: Stopping threshold

    Outputs:
    - Xk: N1 x N2 x K x T solution
    - cgIters: Final number of iterations
    """

    rho1 = params['rho1']
    rho2 = params['rho2']

    # Copy initial guess
    Xk = X_init.copy()

    # Target vector components
    AtB = mat.AtB_ft_Time(A0, B, Bnorms)  # shape: N1 x N2 x K x T
    PtDtZ = mat.PhiTranDiffTran(ZU, Xk.shape)  # shape: N1 x N2 x K x T

    # Initial residual
    Rk = AtB \
        - AtAx(A0, Xk, Bnorms) \
        + rho2 * PtDtZ \
        - rho2 * PtDtDPx(Xk) \
        + rho1 * YV \
        - rho1 * Xk

    Pk = Rk.copy()

    for i in range(params['conjGradIter']):
        Apk = AtAx(A0, Pk, Bnorms) \
            + rho2 * PtDtDPx(Pk) \
            + rho1 * Pk

        RkRk = np.sum(Rk * Rk)
        alphak = RkRk / np.sum(Pk * Apk)

        Xk = Xk + alphak * Pk
        Rkp1 = Rk - alphak * Apk

        if np.linalg.norm(Rkp1) < params['cgEpsilon']:
            break

        betak = np.sum(Rkp1 * Rkp1) / RkRk
        Pk = Rkp1 + betak * Pk
        Rk = Rkp1

    cgIters = i + 1
    return Xk, cgIters

def AtAx(A0ft_stack, X, Bnorms):
    Ax = mat.Ax_ft_Time(A0ft_stack, X, Bnorms)
    return mat.AtB_ft_Time(A0ft_stack, Ax, Bnorms)

def PtDtDPx(X):
    """
    Compute Phi^T D^T D Phi X (temporal difference operator and its adjoint).
    Works for both 1D (N x K x T) and 2D (N1 x N2 x K x T).

    Inputs:
    - X: N x K x T or N1 x N2 x K x T

    Output:
    - PtDtDPx: same shape as X
    """
    dphi = mat.DiffPhiX(X)                  # shape: K x (T-1)
    ptdt = mat.PhiTranDiffTran(dphi, X.shape)
    return ptdt

def convADMM_TGPC(A0, B, X_init, params):
    """
    Generalized ADMM solver for TV-regularized convolutional LASSO using CG.
    Supports both 1D and 2D data.

    Inputs:
    - A0: FFT of dictionary (N x K) or (N1 x N2 x K)
    - B: Data matrix (N x T) or (N1 x N2 x T)
    - X_init: Initial solution (N x K x T) or (N1 x N2 x K x T)
    - params: Dictionary of parameters

    Outputs:
    - X_hat: Estimated sparse coefficients
    - err: Data fidelity error per iteration
    - obj: Objective value per iteration
    - l1_norm: L1 norm per iteration
    - tv_penalty: TV penalty per iteration
    """

    # Unpack parameters
    tolerance = params['tolerance']
    lambda_val = params['lambda']
    gamma = params['gamma']
    rho1 = params['rho1']
    rho2 = params['rho2']
    mu = params['mu']
    adaptRho = params['adaptRho']
    tau = params['tau']
    alpha = params['alpha']
    maxIter = params['maxIter']
    isNonnegative = params['isNonnegative']
    
    is2D = X_init.ndim == 4
    T = X_init.shape[-1]

    # Normalize B and X_init
    Bnorms = np.zeros(T)
    for t in range(T):
        if params['normData']:
            Bnorms[t] = np.linalg.norm(B[..., t])
        else:
            Bnorms[t] = 1
        B[..., t] /= Bnorms[t]
        X_init[..., t] /= Bnorms[t]

    A0 = cp.asarray(A0)
    B = cp.asarray(B)
    Bnorms = cp.asarray(Bnorms)

    # Initialize variables
    Xk = cp.asarray(X_init.copy())
    Xmin = X_init.copy()
    Yk = cp.asarray(X_init.copy())
    Ykp1 = cp.asarray(X_init.copy())
    Vk = cp.asarray(np.zeros_like(Yk))

    Zk = cp.asarray(mat.DiffPhiX(X_init))
    Uk = np.zeros_like(Zk)

    err = np.zeros(maxIter)
    l1_norm = np.zeros(maxIter)
    tv_penalty = np.zeros(maxIter)
    obj = np.zeros(maxIter)

    keep_going = True
    nIter = 0
    count = 0

    while keep_going and (nIter < maxIter):
        nIter += 1

        # CG Solve
        if (nIter > 1) or (np.sum(X_init) == 0):
            if is2D:
                Xkp1, cgIters = conjGrad_TVphi_2D(A0, B, Bnorms, Xk, Yk - Vk, Zk - Uk, params)
            else:
                Xkp1, cgIters = conjGrad_TVphi_1D(A0, B, Bnorms, Xk, Yk - Vk, Zk - Uk, params)
        else:
            Xkp1 = Xk
            cgIters = 0

        # Y-update (soft thresholding)
        for t in range(T):
            Ykp1[..., t] = soft(alpha * Xkp1[..., t] + (1 - alpha) * Yk[..., t] + Vk[..., t], lambda_val[t] / rho1)
        if isNonnegative:
            Ykp1[Ykp1 < 0] = 0
        Vk += alpha * Xkp1 + (1 - alpha) * Yk - Ykp1

    
        # Z-update (TV soft thresholding)
        if gamma > 0:
            Zkp1 = soft(mat.DiffPhiX(Xkp1) + Uk, gamma / rho2)
            Uk += mat.DiffPhiX(Xkp1) - Zkp1

        # Track metrics
        fit = mat.Ax_ft_Time(A0, Xkp1, Bnorms)
        err[nIter - 1] = np.sum((B - fit) ** 2)
        l1_norm[nIter - 1] = sum(lambda_val[t] * np.sum(np.abs(Xkp1[..., t])) for t in range(T))
        tv_penalty[nIter - 1] = gamma * np.sum(np.abs(mat.DiffPhiX(Xkp1)))
        obj[nIter - 1] = 0.5 * err[nIter - 1] + l1_norm[nIter - 1] + tv_penalty[nIter - 1]

        if obj[nIter - 1] <= min(obj):
            Xmin = Xkp1

        if params['verbose']:
            print(f"Iter {nIter}, cgIters {cgIters}, Obj {obj[nIter - 1]:.4e}, "
                  f"Err {0.5 * err[nIter - 1]:.4e}, ||x||_1 {l1_norm[nIter - 1]:.4e}, "
                  f"TVx {tv_penalty[nIter - 1]:.4e}, ||x||_0 {np.sum(Xkp1 > 0)}")

        # Stopping criteria
        if params['stoppingCriterion'] == 'OBJECTIVE_VALUE':
            try:
                keep_going = abs(obj[nIter - 1] - obj[nIter - 2]) > tolerance
            except:
                keep_going = True
        elif params['stoppingCriterion'] == 'COEF_CHANGE':
            diff_x = np.sum(np.abs(Xkp1 - Xk)) / Xk.size
            keep_going = diff_x > tolerance
        else:
            raise ValueError("Undefined stopping criterion.")

        # Adaptive rho update
        if adaptRho:
            skY = rho1 * np.linalg.norm(Ykp1 - Yk)
            rkY = np.linalg.norm(Xkp1 - Ykp1)
            if rkY > mu * skY:
                rho1 *= tau
            elif skY > mu * rkY:
                rho1 /= tau
        
            if gamma > 0:
                skZ = rho2 * np.linalg.norm(Zkp1 - Zk)
                rkZ = np.linalg.norm(mat.DiffPhiX(Xkp1) - Zkp1)
                if rkZ > mu * skZ:
                    rho2 *= tau
                elif skZ > mu * rkZ:
                    rho2 /= tau

        # Early stopping if objective increases
        if (nIter > 10) and (obj[nIter - 2] < obj[nIter - 1]):
            count += 1
            if count > 20:
                keep_going = False
        else:
            count = 0

        # Update iterates
        Xk = Xkp1.copy()
        Yk = Ykp1.copy()
        if gamma > 0:
            Zk = Zkp1.copy()

    X_hat = Xmin.copy()
    if isNonnegative:
        X_hat[X_hat < 0] = 0

    return X_hat, err[:nIter], obj[:nIter], l1_norm[:nIter], tv_penalty[:nIter]


def soft(x, T):
    if np.sum(np.abs(T)) == 0:
        return x
    else:
        y = np.maximum(np.abs(x) - T, 0)
        return np.sign(x) * y

