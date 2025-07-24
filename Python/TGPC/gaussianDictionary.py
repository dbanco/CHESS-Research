# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:45:30 2024

@author: dpqb1
"""
import numpy as np
from scipy.signal import convolve2d

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

def gaussian_basis_wrap_2D(N1, N2, mu1, mu2, sigma, scaling='standard'):
    """
    Create a wrapped 2D Gaussian by outer-product of two 1D wrapped Gaussians.

    Inputs:
    - N1, N2: Dimensions of the 2D grid
    - mu1, mu2: Centers along each dimension (float, not int)
    - sigma: Scalar or (sigma1, sigma2)
    - scaling: '2-norm', '1-norm', 'max', 'rms', or 'standard'

    Output:
    - G: 2D array of shape (N1, N2)
    """
    if np.isscalar(sigma):
        sigma1 = sigma2 = sigma
    else:
        sigma1, sigma2 = sigma

    g1 = gaussian_basis_wrap_1D(N1, mu1, sigma1, scaling='standard')
    g2 = gaussian_basis_wrap_1D(N2, mu2, sigma2, scaling='standard')

    G = np.outer(g1, g2)

    if scaling == '2-norm':
        G /= np.linalg.norm(G)
    elif scaling == '1-norm':
        G /= np.sum(np.abs(G))
    elif scaling == 'max':
        G /= np.max(G)
    elif scaling == 'rms':
        G /= np.sqrt(np.sum(G**2) / (N1 * N2))
    else:  # 'standard': already scaled by scalar Gaussian
        pass

    return G

def peakDictionary(P):
    """
    Unified 1D/2D dictionary builder.
    - 1D if 'N' in P
    - 2D if 'N1' and 'N2' in P with K1, K2, sigmas1, sigmas2
    Returns:
    - D: shape (N, K) for 1D, (N1, N2, K1*K2) for 2D
    """
    if 'N1' in P and 'N2' in P:
        N1, N2 = P['N1'], P['N2']
        K1, K2 = P['K1'], P['K2']
        sigmas1 = P['sigmas1']
        sigmas2 = P['sigmas2']

        D = np.zeros((N1, N2, K1 * K2), dtype=np.float64)
        for k1 in range(K1):
            for k2 in range(K2):
                k = k1 * K2 + k2  # Flattened index
                D[..., k] = gaussian_basis_wrap_2D(
                    N1, N2, mu1=N1 // 2, mu2=N2 // 2,
                    sigma=(sigmas1[k1], sigmas2[k2]), scaling='2-norm'
                )
        return D

    elif 'N' in P:
        N, K = P['N'], P['K']
        sigmas = P['sigmas']
        D = np.zeros((N, K), dtype=np.float64)
        for k in range(K):
            D[:, k] = gaussian_basis_wrap_1D(N, N // 2, sigmas[k], '2-norm')
        return D

    else:
        raise ValueError("peakDictionary: expected keys for 1D ('N') or 2D ('N1','N2','K1','K2').")

def peakDictionaryFFT(P):
    """
    Unified 1D/2D FFT dictionary builder.
    Returns:
    - Dft: shape (N, K) for 1D, (N1, N2, K1*K2) for 2D
    """
    if 'N1' in P and 'N2' in P:
        D = peakDictionary(P)  # shape (N1, N2, K1*K2)
        N1, N2, K = D.shape
        Dft = np.zeros((N1, N2, K), dtype=np.complex128)
        for k in range(K):
            Dft[..., k] = np.fft.fft2(D[..., k])
        return Dft

    elif 'N' in P:
        D = peakDictionary(P)  # shape (N, K)
        return np.fft.fft(D, axis=0)

    else:
        raise ValueError("peakDictionaryFFT: expected keys for 1D ('N') or 2D ('N1','N2','K1','K2').")


def peakDictionarySeparable(P):
    """
    Generate zero-mean Gaussian basis function vectors with unit 2-norm.

    Inputs:
    - P: Dictionary parameters dictionary containing:
        - P.N: signal dimensions
        - P.K: Number of dictionary shapes for each dimensions [K1,K2,...]
        - P.sigmas: Vector of width parameters for each dimension [sigs1,sigs2,...]

    Outputs:
    - D: Dictionary atoms [N, K]
    """

    D = [ np.zeros((P['N'][:],P['K'][k])) for k in len(P['K']) ]
    for nk in len(D):
        for k,sigk in enumerate(P['sigmas'][nk]):
            D[nk][k] = gaussian_basis_wrap_1D(P['N'], np.floor(P['N'] / 2), sigk, '2-norm')

    return D

def peakDictionarySeparableFFT(P):
    """
    Generate zero-mean Gaussian basis function vectors with unit 2-norm.

    Inputs:
    - P: Dictionary parameters dictionary containing:
        - P.N: signal dimensions
        - P.K: Number of dictionary shapes for each dimensions [K1,K2,...]
        - P.sigmas: Vector of width parameters for each dimension [sigs1,sigs2,...]

    Outputs:
    - D: Dictionary atoms [N, K]
    """

    D = [ np.zeros((P['N'][:],P['K'][k])) for k in len(P['K']) ]
    for nk in len(D):
        for k,sigk in enumerate(P['sigmas'][nk]):
            D[nk][k] = np.fft.fft(gaussian_basis_wrap_1D(P['N'], np.floor(P['N'] / 2), sigk, '2-norm'))

    return D

def computeAWMV(x, sigmas1, sigmas2=None):
    """
    Compute amplitude-weighted mean variance (AWMV) for 1D or 2D data.

    Parameters:
    - x: ndarray
        Shape [N, K, T] for 1D
        Shape [N1, N2, K1*K2, T] for 2D
    - sigmas1: array of shape (K,) for 1D or (K1,) for 2D
    - sigmas2: None for 1D or (K2,) for 2D

    Returns:
    - awmv: float for 1D, or (2,) array [awmv_x, awmv_y] for 2D
    """
    x = np.array(x)

    if sigmas2 is None:
        # 1D case
        sigmas1 = np.array(sigmas1).reshape(1, -1)  # shape (1, K)
        amp = np.sum(x, axis=tuple(range(x.ndim - 2)))  # sum over N, T -> shape (K,)
        total = np.sum(amp)

        if total == 0:
            return 0.0

        awmv = np.sum(sigmas1 * amp.reshape(1, -1) / total)
        return awmv

    else:
        # 2D case
        K1, K2 = len(sigmas1), len(sigmas2)
        sigmas1 = np.array(sigmas1).reshape(1, K1)  # shape (1, K1)
        sigmas2 = np.array(sigmas2).reshape(1, K2)  # shape (1, K2)

        # Reshape x from [N1, N2, K1*K2, T] to [N1, N2, K1, K2, T]
        x = x.reshape(*x.shape[:2], K1, K2, -1)

        # Sum over N1, N2, T → axes 0, 1, 4
        amp1 = np.sum(x, axis=(0, 1, 3))  # shape (K1,T)
        amp2 = np.sum(x, axis=(0, 1, 2))  # shape (K2,T)
        total = np.sum(amp1)

        if total == 0:
            return np.array([0.0, 0.0])

        awmv1 = np.dot(sigmas1,amp1) / total  # weighted sum over K1
        awmv2 = np.dot(sigmas2,amp2) / total  # weighted sum over K2

        return np.array([awmv1, awmv2])

def generateExampleData(N, T):
    """
    Generate example 1D Gaussian peaks with time-varying mean and width,
    and Poisson noise.
    """
    numSpots = 2
    B = np.zeros((N, T))
    B_noise = np.zeros((N, T))
    amplitude = 80 * np.array([0.4, 0.7]) + 1

    t_vals = np.linspace(0, 1, T)
    # Evolving means, staying inside the domain
    mean_base = N * np.array([0.3, 0.7])
    mean_amp = N * 0.1
    mean_param = np.array([mean_base[0] + mean_amp * t_vals**2,
                           mean_base[1] + mean_amp * 2*t_vals**2])

    # Evolving widths
    widths = np.array([1.5 + (3*t_vals)**2,
                       1 + 0.5 * (2*t_vals)**2])

    awmv_true = np.zeros(T)
    for t in range(T):
        for i in range(numSpots):
            b = gaussian_basis_wrap_1D(N,
                                       mean_param[i, t],
                                       widths[i, t],
                                       '2-norm')
            awmv_true[t] += amplitude[i] * widths[i, t]
            B[:, t] += amplitude[i] * b

        awmv_true[t] /= np.sum(amplitude)
        B_noise[:, t] = np.random.poisson(B[:, t])

    return B, B_noise, awmv_true


def generateExampleData2(dims):
    """
    Generate 2D Gaussian peaks with time-varying mean and width, and Poisson noise.
    """
    numSpots = 2
    N1, N2, T = dims
    B = np.zeros(dims)
    B_noise = np.zeros(dims)
    amplitude = 80 * np.array([0.4, 0.7]) + 1

    t_vals = np.linspace(0, 1, T)

    # Evolving means (bounded inside image)
    mean_base = np.array([[N1 * 0.3, N1 * 0.6],   # Row direction
                          [N2 * 0.3, N2 * 0.7]])  # Column direction

    mean_amp = np.array([[N1 * 0.1, N1 * 0.05],   # Row direction
                         [N2 * 0.1, N2 * 0.05]])  # Column direction

    mean_param = np.zeros((2, numSpots, T))  # [axis, spot, time]
    mean_param[0, 0, :] = mean_base[0, 0] + mean_amp[0, 0] * t_vals**2
    mean_param[0, 1, :] = mean_base[0, 1] + mean_amp[0, 1] * 2 * t_vals**2
    mean_param[1, 0, :] = mean_base[1, 0] + mean_amp[1, 0] * t_vals**2
    mean_param[1, 1, :] = mean_base[1, 1] + mean_amp[1, 1] * 2 * t_vals**2

    # Evolving widths
    widths = np.zeros((2, numSpots, T))  # [axis, spot, time]
    widths[0, 0, :] = 1.5 + (3 * t_vals)**2
    widths[0, 1, :] = 1.0 + 0.5 * (2 * t_vals)**2
    widths[1, 0, :] = 2.0 + (2.5 * t_vals)**2
    widths[1, 1, :] = 1.2 + 0.4 * (2 * t_vals)**2

    awmv_true = np.zeros((T, 2))  # [time, axis]
    for t in range(T):
        for i in range(numSpots):
            b1 = gaussian_basis_wrap_1D(N1,
                                        mean_param[0, i, t],
                                        widths[0, i, t],
                                        '2-norm')
            b2 = gaussian_basis_wrap_1D(N2,
                                        mean_param[1, i, t],
                                        widths[1, i, t],
                                        '2-norm')
            outer = amplitude[i] * convolve2d(b1[:, None], b2[None, :])
            B[:, :, t] += outer

            awmv_true[t, 0] += amplitude[i] * widths[0, i, t]
            awmv_true[t, 1] += amplitude[i] * widths[1, i, t]

        awmv_true[t, :] /= np.sum(amplitude)
        B_noise[:, :, t] = np.random.poisson(B[:, :, t])

    return B, B_noise, awmv_true

