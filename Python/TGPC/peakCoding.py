# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:48:32 2024

@author: dpqb1
"""

import numpy as np
import matrixOpsTGPC as mat

# def auto_noise_estimation(B,)

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
    AtB = mat.AtB_ft_1D_Time(A0, B, Bnorms)
    PtDtZ = mat.PhiTranDiffTran_1D(ZU, N)

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


def AtAx(A0ft_stack, X, Bnorms):
    return mat.AtB_ft_1D_Time(A0ft_stack, mat.Ax_ft_1D_Time(A0ft_stack, X, Bnorms), Bnorms)


def PtDtDPx(X):
    N = X.shape[0]
    return mat.PhiTranDiffTran_1D(mat.DiffPhiX_1D(X), N)


def convADMM_LASSO_CG_TVphi_1D(A0, B, X_init, params):
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

    N, K, T = X_init.shape
    Bnorms = np.zeros(T)
    for j in range(T):
        if params['normData']:
            Bnorms[j] = np.linalg.norm(B[:, j])
        else:
            Bnorms[j] = 1
        B[:, j] /= Bnorms[j]
        X_init[:, :, j] /= Bnorms[j]

    Xk = X_init.copy()
    Xmin = X_init.copy()

    Yk = X_init.copy()
    Ykp1 = X_init.copy()
    Vk = np.zeros_like(Yk)

    Zk = mat.DiffPhiX_1D(X_init)
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

        if (nIter > 1) or (np.sum(X_init) == 0):
            Xkp1, cgIters = conjGrad_TVphi_1D(A0, B, Bnorms, Xk, (Yk - Vk), (Zk - Uk), params)
        else:
            Xkp1 = Xk
            cgIters = 0

        for t in range(T):
            Ykp1[:, :, t] = soft(alpha * Xkp1[:, :, t] + (1 - alpha) * Yk[:, :, t] + Vk[:, :, t], lambda_val[t] / rho1)
        if isNonnegative:
            Ykp1[Ykp1 < 0] = 0
        Vk += alpha * Xkp1 + (1 - alpha) * Yk - Ykp1

        Zkp1 = soft(mat.DiffPhiX_1D(Xkp1) + Uk, gamma / rho2)
        Uk += mat.DiffPhiX_1D(Xkp1) - Zkp1

        fit = mat.Ax_ft_1D_Time(A0, Xkp1, Bnorms)
        err[nIter - 1] = np.sum((B - fit) ** 2)
        Xsum = 0
        for t in range(T):
            Xsum += lambda_val[t] * np.sum(np.abs(Xkp1[:, :, t]))
        l1_norm[nIter - 1] = Xsum
        tv_penalty[nIter - 1] = gamma * np.sum(np.abs(mat.DiffPhiX_1D(Xkp1)))
        obj[nIter - 1] = 0.5 * err[nIter - 1] + l1_norm[nIter - 1] + tv_penalty[nIter - 1]

        if obj[nIter - 1] <= min(obj):
            Xmin = Xkp1

        if params['verbose']:
            print('Iter', nIter, 'cgIters', cgIters, 'Rho1', rho1, 'Rho2', rho2, 'Obj', obj[nIter - 1],
                  'Err', 0.5 * err[nIter - 1], '||x||_1', l1_norm[nIter - 1], 'TVx', tv_penalty[nIter - 1],
                  '||x||_0', np.sum(Xkp1 > 0))

        if params['plotProgress']:
            # Plot progress if needed
            pass

        # Check stopping criterion
        if params['stoppingCriterion'] == 'OBJECTIVE_VALUE':
            try:
                criterionObjective = abs(obj[nIter - 1] - obj[nIter - 2])
                keep_going = (criterionObjective > tolerance)
            except:
                keep_going = True
        elif params['stoppingCriterion'] == 'COEF_CHANGE':
            diff_x = np.sum(np.abs(Xkp1 - Xk)) / Xk.size
            keep_going = (diff_x > tolerance)
        else:
            raise ValueError('Undefined stopping criterion.')

        if adaptRho:
            skY = rho1 * np.linalg.norm(Ykp1 - Yk)
            skZ = rho2 * np.linalg.norm(Zkp1 - Zk)
            rkY = np.linalg.norm(Xkp1 - Ykp1)
            diff_rkZ = mat.DiffPhiX_1D(Xkp1) - Zkp1
            rkZ = np.linalg.norm(diff_rkZ)
            if rkY > mu * skY:
                rho1 *= tau
            elif skY > mu * rkY:
                rho1 /= tau
            if rkZ > mu * skZ:
                rho2 *= tau
            elif skZ > mu * rkZ:
                rho2 /= tau

        if (nIter > 10) and (obj[nIter - 2] < obj[nIter - 1]):
            count += 1
            if count > 20:
                keep_going = False
        else:
            count = 0

        Xk = Xkp1.copy()
        Yk = Ykp1.copy()
        Zk = Zkp1.copy()

    X_hat = Xmin.copy()
    if isNonnegative:
        X_hat[X_hat < 0] = 0

    err = err[:nIter]
    obj = obj[:nIter]

    return X_hat, err, obj, l1_norm, tv_penalty

def soft(x, T):
    if np.sum(np.abs(T)) == 0:
        return x
    else:
        y = np.maximum(np.abs(x) - T, 0)
        return np.sign(x) * y

