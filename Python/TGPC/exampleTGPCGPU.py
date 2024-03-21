# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:55:22 2024

@author: dpqb1
"""

import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
import matrixOpsTGPCGPU as mat 
import gaussianDictionary as gd
import peakCodingGPU as pc

# Generate example data
T = 10
N = 101
K = 20;
B, B_poiss, awmv_true = gd.generateExampleData(N, T)

# Define parameters
P = {}
P['N'] = B.shape[0]
P['K'] = K
P['sigmas'] = np.linspace(1/2, 16, P['K'])
params = {
    'adaptRho': 1,
    'mu': 2,
    'tau': 1.05,
    'alpha': 1.8,
    'isNonnegative': 1,
    'normData': 1,
    'stoppingCriterion': 'OBJECTIVE_VALUE',
    'maxIter': 100,
    'conjGradIter': 100,
    'tolerance': 1e-8,
    'cgEpsilon': 1e-6,
    'plotProgress': 0,
    'verbose': 1
}
P['params'] = params

# Construct dictionary
A0ft = cp.asarray(gd.peakDictionaryFFT(P))
A0 = gd.peakDictionary(P)

# Setup and solve
params['rho1'] = 1
params['lambda'] = 4e-4 * np.ones(T)
params['rho2'] = 0.1
params['gamma'] = 0
admmOutIndep = pc.convADMM_LASSO_CG_TVphi_1D(A0ft, B, np.zeros((N, P['K'], T)), params)
X_hat_indep = admmOutIndep[0];
B_hat_indep = mat.Ax_ft_1D_Time(A0ft, X_hat_indep)

params['gamma'] = 5e-3
admmOut = pc.convADMM_LASSO_CG_TVphi_1D(A0ft, B_poiss, np.zeros((N, P['K'], T)), params)
X_hat = admmOut[0]
B_hat = mat.Ax_ft_1D_Time(A0ft, X_hat)

# Plot awmv recovery
fig, axs = plt.subplots(1, 4, figsize=(15, 4))

# Plot truth
axs[0].imshow(B.T, aspect='auto', origin='lower', cmap='viridis')
axs[0].set_title('Truth')

# Plot independent reconstruction
axs[1].imshow(B_hat_indep.T, aspect='auto', origin='lower', cmap='viridis')
err1a = np.linalg.norm(B_hat_indep - B) / np.linalg.norm(B)
err1b = np.linalg.norm(B_hat_indep - B_poiss) / np.linalg.norm(B_poiss)
axs[1].set_title(f'Recon (Indep), err: {err1a:.3f}, {err1b:.3f}')

# Plot coupled reconstruction
axs[2].imshow(B_hat.T, aspect='auto', origin='lower', cmap='viridis')
err2a = np.linalg.norm(B_hat - B) / np.linalg.norm(B)
err2b = np.linalg.norm(B_hat - B_poiss) / np.linalg.norm(B_poiss)
axs[2].set_title(f'Recon (Coupled), err: {err2a:.3f}, {err2b:.3f}')

# Plot Poisson data
axs[3].imshow(B_poiss.T, aspect='auto', origin='lower', cmap='viridis')
axs[3].set_title('Poisson Data')

fig4, ax4 = plt.subplots(figsize=(8, 6))
awmv1 = gd.computeAWMV_1D(X_hat_indep, P['sigmas'])
awmv2 = gd.computeAWMV_1D(X_hat, P['sigmas'])
awmv_err3 = np.linalg.norm(awmv1 - awmv_true) / np.linalg.norm(awmv_true)
awmv_err4 = np.linalg.norm(awmv2 - awmv_true) / np.linalg.norm(awmv_true)
ax4.plot(awmv_true, label='True', linewidth=2)
ax4.plot(awmv1, label=f'Indep: {awmv_err3:.3f}', linewidth=2)
ax4.plot(awmv2, label=f'Coupled: {awmv_err4:.3f}', linewidth=2)
ax4.set_ylabel('AWMV')
ax4.set_xlabel('t')
ax4.set_title('AWMV')
ax4.legend(loc='best')
plt.show()
