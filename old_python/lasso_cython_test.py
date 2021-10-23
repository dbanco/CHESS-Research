# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:56:01 2016

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
stds = [5, 10, 7]

np.random.normal(0,1,num_theta)

noise = np.random.normal(0,1,num_theta)
signal = noise

for i in range(3):
    signal += scale[i]*EM.gaussian_basis_wrap(num_theta,
                                           dtheta,
                                           means[i],
                                           (stds[i]*dtheta)**2)

#%% Create Dictionary
num_var = 15

# Local maxima
maxima = argrelmax(signal)[0]

# Variances
test_vars = np.linspace((dtheta),(np.pi/16),num_var)**2
variances = test_vars

B = np.empty([num_theta,0])
B0 = np.empty([num_theta,0])

for i, var in enumerate(variances):  
    basis_vec = EM.gaussian_basis_wrap(num_theta,dtheta,0,var)
    B0 = np.hstack([B0,basis_vec[:,None]])
    
for shift in maxima:
    B_shift = np.vstack([B0[-shift:,:], B0[0:-shift,:]])
    B = np.hstack([B,B_shift])


B = np.array(B, order='F', copy=True)

#%% Test Lasso fitting routines
import lasso_coord_ascent as coord
reload(coord)
alpha = 0.08
max_iters = 1000

a = time.time()
coefs1 = coord.lasso_solver(B,signal,alpha,max_iters,positive=1)
y_hat1 = B.dot(coefs1)
b = time.time()
cy_time = b-a

a = time.time()
#coefs2, y_hat2, rel_error = RingIP.lasso_solver(B,signal,l1_ratio,max_iters)
b = time.time()
py_time = b-a

a = time.time()
lasso = Lasso(alpha=alpha,
              max_iter=max_iters,
              fit_intercept=0,
              positive=True)
lasso.fit(B,signal)
y_hat3 = B.dot(lasso.coef_)
b = time.time()  
sci_time = b-a

print('py_time')
print(py_time)
print('cy_time')
print(cy_time)
print('sci_time')
print(sci_time)



#%%
#norm(B.dot(lasso.coef_) - signal)/norm(signal)
#norm(B.dot(coefs1) - signal)/norm(signal)
#norm(B.dot(coefs2) - signal)/norm(signal)

plt.plot(signal,'-og')
plt.plot(y_hat1,'-ob')
#plt.plot(y_hat2,'-om')
plt.plot(y_hat3,'-or')