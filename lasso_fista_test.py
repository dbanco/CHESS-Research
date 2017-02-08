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
from math import sqrt

import RingImageProcessing as RingIP
import EllipticModels as EM
reload(RingIP)
reload(EM)


def soft_thresh(x, l):
    return np.sign(x) * np.maximum(np.abs(x) - l, 0.)    
    
def fista(A, b, l, maxit):
    x = np.zeros(A.shape[1])
    t = 1
    z = x.copy()
    L = norm(A)**2
    for it in range(maxit):
        xold = x.copy()
        z = z + A.T.dot(b - A.dot(z))/L
        x = soft_thresh(z,1/L)
        t0 = t
        t = (1. + sqrt(1. + 4. * t ** 2)) / 2.
        z = x + ((t0 - 1.) / t) * (x - xold)

    return x

def Ax_ft(A0ft, x_ft):  
    # Inputs:
    #    A0ft    Each column is the first column of a circulant matrix Ai
    #    x_ft    Each column is the fourier transform of the coefficients
    #            corresponding to the circulant matrix Ai
    Ax = np.zeros((A0ft.shape[0]))
    
    for ii in range(A0ft.shape[1]):
        Ax += np.fft.ifft(np.multiply(A0ft[:,ii],x_ft[:,ii])).real
        
    return Ax
    
def ATb_ft(A0_tran_ft,R,num_var,num_maxima):
    # Inputs:
    #    A0_tran_ft    Each column is the first row of a circulant matrix Ai
    #    R             Residual vector 
    # Output:
    #    ATb           Remember that this will be zero padded

    ATb = np.zeros((num_var*num_maxima))
    R_ft = np.fft.fft(R)        
    
    print(A0_tran_ft.shape)    
    
    # num_vars    
    for ii in range(num_var):
        ATb[ii*num_maxima:(ii+1)*num_maxima] = np.fft.ifft(np.multiply(A0_tran_ft[:,ii],R_ft)).real[maxima]   

    return ATb
    
def x_to_x_ft(x, maxima, m, n):
    # x coefficient vector
    # m number of rows in x_ft
    # n number of columns in x_ft


    x_pad = np.zeros((m,n))
    x_ft = np.zeros((m,n))
    num_maxima = len(maxima)
    for i in range(n):
        x_pad[maxima,i] = x[i*num_maxima:(i+1)*num_maxima]      
        
    x_ft = np.fft.fft(x_pad,axis=0)
    
    return x_ft
    
def fista_circulant(A, A0, A0_tran, maxima, b, l, maxit):
    m = A.shape[0]    
    n = A0.shape[1]
    A0ft = np.fft.fft(A0,axis=0)    
    A0_tran_ft = np.fft.fft(A0_tran,axis=0)  
    
    x = np.zeros(A.shape[1])
    t = 1
    z = x.copy()
    L = norm(A)**2
    for it in range(maxit):
        xold = x.copy()
        
        # Arrange x coefficents as matrix in fourier domain 
        z_ft = x_to_x_ft(z, maxima, m, n)
        R = b - Ax_ft(A0ft, z_ft)
        z = z + ATb_ft(A0_tran_ft,R,n,len(maxima))/L
        #z = z + A.T.dot(R)/L
        x = soft_thresh(z,1/L)
        t0 = t
        t = (1. + sqrt(1. + 4. * t ** 2)) / 2.
        z = x + ((t0 - 1.) / t) * (x - xold)

    return x
    
    
#%% Generate Data

num_theta = 250
dtheta = 2*np.pi/num_theta

scale = [200, 400, 330]
means = [40, 140, 200]
stds = [5, 10, 7]

np.random.normal(0,1,num_theta)

noise = np.random.normal(0,1,num_theta)
signal = np.zeros(num_theta)

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
B0_tran = np.zeros([num_theta,num_var])

for i, var in enumerate(variances):  
    basis_vec = EM.gaussian_basis_wrap(num_theta,dtheta,0,var)
    B0 = np.hstack([B0,basis_vec[:,None]])
    for shift in maxima:
        B_shift = np.concatenate([basis_vec[-shift:], basis_vec[0:-shift]])
        B = np.hstack([B,B_shift[:,None]]) 
    
B_tran = np.transpose(B)

for i, var in enumerate(variances): 
    B0_tran[0,i] = B0[0,i]
    B0_tran[1::,i] = np.flipud(B0[1::,i])

B0ft = np.fft.fft(B0,axis=0)
B0_tranft = np.fft.fft(B0_tran,axis=0)

B = np.array(B, order='F', copy=True)
B_tran = np.array(B_tran, order='F', copy=True)
B0ft = np.array(B0ft, order='F',copy=True)
B0_tranft = np.array(B0_tranft, order='F',copy=True)
B0 = np.array(B0, order='F',copy=True)
B0_tran = np.array(B0_tran, order='F',copy=True)

#%% Test Lasso fitting routines
import circulant_lasso_ca_np as coord
reload(coord)
alpha = 0.08
max_iters = 50

#a = time.time()
#coefs1 = coord.lasso_solver(B,B0,B0_tran,maxima,signal,alpha,max_iters,positive=1)
#y_hat1 = B.dot(coefs1)
#b = time.time()
#cy_time = b-a

a = time.time()
coefs2 = lasso_solver(B,B_tran,B0,B0_tran,maxima,signal,alpha,max_iters,positive=1)
y_hat2 = B.dot(coefs2)
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

a = time.time()
coefs4 = lasso_solver_batch(B,B_tran,B0,B0_tran,maxima,signal,alpha,max_iters,positive=1)
y_hat4 = B.dot(coefs4)
b = time.time()
batch_time = b-a

a = time.time()
coefs5 = fista(B,signal,alpha,max_iters)
y_hat5 = B.dot(coefs5)
b = time.time()
fista_time = b-a

a = time.time()
coefs6 = fista_circulant(B,B0,B0_tran,maxima,signal,alpha,max_iters)
y_hat6 = B.dot(coefs6)
b = time.time()
circ_time = b-a

print('py_time')
print(py_time)
#print('cy_time')
#print(cy_time)
print('sci_time')
print(sci_time)
print('batch_time')
print(batch_time)
print('fista_time')
print(fista_time)
print('circ_time')
print(circ_time)
#%%
#norm(B.dot(lasso.coef_) - signal)/norm(signal)
#norm(B.dot(coefs1) - signal)/norm(signal)
#norm(B.dot(coefs2) - signal)/norm(signal)

plt.plot(signal,'-og')
#plt.plot(y_hat1,'-ob')
#plt.plot(y_hat2,'-xm')
plt.plot(y_hat3,'-or')
#plt.plot(y_hat4,'-sk')
plt.plot(y_hat5,'-sb')
plt.plot(y_hat6,'-xm')


