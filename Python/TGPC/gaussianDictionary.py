# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:45:30 2024

@author: dpqb1
"""
import numpy as np


def gaussian_basis_wrap_1D(N, mu, sigma, scaling='standard'):
    """
    Generate a Gaussian peak function vector.

    Inputs:
    - N: Vector length
    - mu: Mean of the Gaussian basis function
    - sigma: Standard deviation of the Gaussian basis function
    - scaling: Scaling type (default is 'standard')
               Options:
                   '2-norm': Unit 2-norm scaling
                   '1-norm': Unit 1-norm scaling
                   'max': Unit max scaling factor
                   'rms': Unit root-mean-square

    Outputs:
    - b: N x 1 vector
    """

    # Compute theta distances with wrapping at boundaries
    idx = np.arange(1, N + 1)
    wrapN = lambda x, N: int(1 + np.mod(x - 1, N))
    opposite = (idx[wrapN(np.floor(mu - N / 2), N) - 1] +
                idx[wrapN(np.ceil(mu - N / 2), N) - 1]) / 2
    if opposite == mu:
        opposite = 0.5
    dist1 = np.abs(mu - idx)
    dist2 = N / 2 - np.abs(opposite - idx)
    dist = np.minimum(dist1, dist2)
    dist_sq_theta = dist**2  # num_theta length vector

    # Compute values
    b = np.exp(-dist_sq_theta / (2 * sigma**2))
    if scaling == '2-norm':
        b /= np.linalg.norm(b)
    elif scaling == '1-norm':
        b /= np.sum(np.abs(b))
    elif scaling == 'max':
        b /= np.max(b)
    elif scaling == 'rms':
        b /= np.sqrt(np.sum(b**2) / N)
    else:
        b /= (sigma * np.sqrt(2 * np.pi))

    return b

def peakDictionary(P):
    """
    Generate zero-mean Gaussian basis function vectors with unit 2-norm.

    Inputs:
    - P: Dictionary parameters dictionary containing:
        - P.N: 1D signal length
        - P.K: Number of dictionary entries
        - P.sigmas: Vector of width parameters for dictionary

    Outputs:
    - D: Dictionary atoms [N, K]
    """

    D = np.zeros((P['N'], P['K']))
    for k in range(P['K']):
        D[:, k] = gaussian_basis_wrap_1D(P['N'], np.floor(P['N'] / 2), P['sigmas'][k], '2-norm')

    return D

def peakDictionaryFFT(P):
    """
    Generate FFT of zero mean Gaussian basis function vectors with unit 2-norm.

    Inputs:
    - P: Dictionary parameters dictionary containing:
        - P.N: 1D signal length
        - P.K: Number of dictionary entries
        - P.sigmas: Vector of width parameters for dictionary

    Outputs:
    - D: FFT of dictionary atoms [N, K]
    """

    D = np.zeros((P['N'], P['K']), dtype=np.complex128)
    for k in range(P['K']):
        D[:, k] = np.fft.fft(gaussian_basis_wrap_1D(P['N'], np.floor(P['N'] / 2), P['sigmas'][k], '2-norm'))

    return D

def computeAWMV_1D(x, sigmas):
    # Inputs:
    # x - (N x K) array of fitted coefficients
    # sigmas - (K x 1) array of dictionary width parameters
    
    # Outputs:
    # awmv - amplitude weighted mean variance (azimuthal)
    
    K = len(sigmas)
    
    # Ensure sigmas is a row vector
    sigmas = np.array(sigmas).reshape(1, -1)
    
    sigma_signal = np.sum(x, axis=0)
    sigma_signal = sigma_signal.T if sigma_signal.shape[0] == K else sigma_signal
    total = np.sum(sigma_signal, axis=0)
    
    awmv = np.sum(sigmas * sigma_signal / total, axis=1)
    
    return awmv

def generateExampleData(N, T):
    """
    Generate example Poisson measurements of Gaussian intensity peaks.
    """
    numSpots = 2
    B = np.zeros((N, T))
    B_noise = np.zeros((N, T))
    amplitude = 80 * np.array([0.4, 0.7]) + 1
    mean_param = N * np.array([0.3, 0.7])
    widths = np.array([5, 8])

    awmv_true = np.zeros(T)
    for t in range(T):
        for i in range(numSpots):
            b = gaussian_basis_wrap_1D(N,
                                        mean_param[i],
                                        widths[i],
                                        '2-norm')
            awmv_true[t] += amplitude[i] * widths[i]
            B[:, t] += amplitude[i] * b

        awmv_true[t] /= np.sum(amplitude)
        B_noise[:, t] = np.random.poisson(B[:, t])

    return B, B_noise, awmv_true