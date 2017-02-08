# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 12:03:01 2016

@author: Dan
"""

## Init
from __future__ import division

import numpy as np
import os
import time

from sklearn.linear_model import Lasso
from numpy.linalg import norm 
import matplotlib.pyplot as plt
from scipy.signal import argrelmax

import RingImageProcessing as RingIP
import EllipticModels as EM
reload(RingIP)
reload(EM)

#%% Generate Data

num_theta = 250
dtheta = 2*np.pi/num_theta

scale = [200, 400, 330]
means = [40, 140, 200]


np.random.normal(0,1,num_theta)

noise = np.random.normal(0,1,num_theta)
signal = 0

i = 1
signal += EM.gaussian_basis_wrap(num_theta,dtheta,0,(10*dtheta)**2)

coefs = np.zeros((num_theta))

for i in range(3):
    coefs[means[i]] = scale[i]



import cython_test as ct

a_fft = np.fft.fft(signal)    
b_fft = np.fft.fft(coefs)

# Convolution via cython and blas
t = time.time()
for l in range(1000):
    a_fft_re = np.real(a_fft)
    a_fft_im = np.imag(a_fft)
    b_fft_re = np.real(b_fft)
    b_fft_im = np.imag(b_fft)
    
    c_fft_re, c_fft_im = ct.element_multiply(a_fft_re, a_fft_im, b_fft_re, b_fft_im)
    
    c_fft = np.zeros((num_theta),dtype=np.complex128)
    c_fft.real = c_fft_re
    c_fft.imag = c_fft_im
    
    c1 = np.fft.ifft(c_fft)
t1 = time.time() - t

# Convolution in python
t = time.time()
for l in range(1000):
    c2_fft = np.multiply(a_fft,b_fft)
    c2 = np.fft.ifft(c2_fft)`
t2 = time.time() - t

# Convolution with re and im separate
t = time.time()
for l in range(1000):
    c3_fft = np.zeros((num_theta),dtype=np.complex128)
    
    c3_fft_re = np.multiply(a_fft.real,b_fft.real)
    c3_fft_re -= np.multiply(a_fft.imag,b_fft.imag)
    c3_fft_im = np.multiply(a_fft.real,b_fft.imag)
    c3_fft_im += np.multiply(a_fft.imag,b_fft.real)
    
    c3_fft = np.zeros((num_theta),dtype=np.complex128)
    c3_fft.real = c3_fft_re
    c3_fft.imag = c3_fft_im
    
    c3 = np.fft.ifft(c3_fft)
t3 = time.time() - t

plt.figure(1)
plt.plot(signal,'o-r')

plt.figure(2)
plt.plot(c1,'s-r',markersize=9)
plt.plot(c2,'o-b')
plt.plot(c3,'x-g')

print([t1,t2,t3])

