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
    
def Ax_ft(Aft, x_ft):
    
    #   Aft
    #   Each column corresponds to a circulant matrix Ai
    #   x_ft
    #   Each column of x_ft is the fourier transform of the coefficients
    #   corresponding to the circulant matrix Ai
    Ax = np.zeros((Aft.shape[0]))
    
    for ii in range(Aft.shape[1]):
        Ax += np.fft.ifft(np.multiply(Aft[:,ii],x_ft[:,ii])).real
        
    return Ax

def x_to_x_ft(x, maxima, m, n):
    # x coefficient vector
    # m number of rows in x_ft
    # n number of columns in x_ft


    x_pad = np.zeros((m,n))
    x_ft = np.zeros((m,n))
    num_maxima = len(maxima)
    for i in range(n):
        x_pad[maxima,i] = x[i*num_maxima:(i+1)*num_maxima]      
    
    plt.close(50)
    plt.figure(50)
    plt.imshow(x_pad,interpolation='nearest')
    plt.colorbar()

    plt.close(51)
    plt.figure(51)
    plt.plot(np.absolute(np.fft.fft(x_pad[:,0])))    
    
    x_ft = np.fft.fft(x_pad,axis=0)
    
    return x_ft
    
    
def fista_compare(A, A0, maxima, b, l, maxit):
    m = A.shape[0]    
    n = A0.shape[1]
    Aft = np.fft.fft(A0,axis=0)    
    
    x1 = np.zeros(A.shape[1])
    x2 = np.zeros(A.shape[1])
    
    t = 1
    z1 = x1.copy()
    z2 = x2.copy()
    L = norm(A)**2
    for it in range(maxit):
        xold1 = x1.copy()
        xold2 = x2.copy()

        plt.figure(it)
        plt.plot(z1,'bo-')        
        plt.plot(z2,'rs-')       
      
        z1 = z1 + A.T.dot(b - A.dot(z1))/L
        
        # Arrange x coefficents as matrix in fourier domain 
        z_ft = x_to_x_ft(z2, maxima, m, n)
        z2 = z2 + A.T.dot(b - Ax_ft(Aft, z_ft))/L
        
        x1 = soft_thresh(z1,1/L)
        x2 = soft_thresh(z2,1/L)
        
        plt.close(30)
        plt.figure(30)
        plt.plot(Ax_ft(Aft, x_to_x_ft(x2,maxima,m,n)),'gs-')
        
        t0 = t
        t = (1. + sqrt(1. + 4. * t ** 2)) / 2.
        z1 = x1 + ((t0 - 1.) / t) * (x1 - xold1)
        z2 = x2 + ((t0 - 1.) / t) * (x2 - xold2)
        print(norm(z1-z2))  
    return x1,x2

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
B0_tran = np.zeros([num_var,num_theta])

for i, var in enumerate(variances):  
    basis_vec = EM.gaussian_basis_wrap(num_theta,dtheta,0,var)
    B0 = np.hstack([B0,basis_vec[:,None]])
    for shift in maxima:
        B_shift = np.concatenate([basis_vec[-shift:], basis_vec[0:-shift]])
        B = np.hstack([B,B_shift[:,None]]) 
    
B_tran = np.transpose(B)

for i, var in enumerate(variances): 
    B0_tran[i,maxima] = B[0,i*len(maxima):(i+1)*len(maxima)]

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
max_iters = 3

coefs1, coefs2 = fista_compare(B,B0,maxima,signal,alpha,max_iters)

y_hat1 = B.dot(coefs1)
y_hat2 = B.dot(coefs2)

coefs_ft = x_to_x_ft(coefs1,maxima,num_theta,num_var)

#plt.figure(51)
#plt.imshow(np.absolute(coefs_ft),interpolation='nearest')
#plt.colorbar()

y_hat_test = Ax_ft(B0ft,coefs_ft)
#%%
#norm(B.dot(lasso.coef_) - signal)/norm(signal)
#norm(B.dot(coefs1) - signal)/norm(signal)
#norm(B.dot(coefs2) - signal)/norm(signal)

plt.figure(1000)
plt.plot(signal,'-og')
plt.plot(y_hat1,'-sb')
plt.plot(y_hat2,'-om')
plt.plot(y_hat_test,'-xr')

plt.figure(20)
plt.imshow(B,interpolation='nearest')


plt.figure(21)
plt.imshow(B0,interpolation='nearest')