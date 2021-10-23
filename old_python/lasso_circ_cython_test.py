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


def lasso_solver_batch(A,
                       A_tran,
                 A0,
                 A0_tran,
                 maxima,
                 b,
                 alpha, 
                 max_iter,
                 positive=0):

    m = A.shape[0]
    n = A.shape[1]
    num_vars = A0.shape[1]
    num_maxima = maxima.shape[0]

    # precompute sum of sqaures of columns
    z = np.sum(A**2,0) 
   
    # initialize x
    x = np.zeros(n)
    x_old = np.zeros(num_maxima)             
    AtR = np.zeros(m)
    R = np.zeros(m) 


    rho = np.zeros(num_maxima)   
    tol = 1e-4
    gap = tol + 1.0
    d_x_tol = tol       
       
    # R = b - np.dot(A,x) but x is zeros   
    R = b.copy()
    
    # tol *= np.dot(b,b)
    tol *= np.dot(b,b)
 
    for iters in range(max_iter):
        x_max = 0.0
        d_x_max = 0.0
        
        for ii in range(num_vars):
            # 1)            
            x_old = x[ii*num_maxima:(ii+1)*num_maxima]
            
            # Remove contribution to residual of current coefficient
            # R += x_ii*A[:,ii]
            R += np.dot(A[:,ii*num_maxima:(ii+1)*num_maxima],x[ii*num_maxima:(ii+1)*num_maxima])
#            plt.figure(1)
#            plt.plot(R)
#            plt.pause(0.25)
            # 2)
            # Soft threshold update with optional positivity constraint
            rho = np.dot(A_tran[ii*num_maxima:(ii+1)*num_maxima,:],R)        
#            plt.figure(2)
#            plt.plot(rho,'o')
#            plt.pause(0.25)            
            # rho = np.dot(A[:,ii],R)
            if positive:
                rho[rho<0]=0
            x[num_maxima*ii:num_maxima*(ii+1)] = np.sign(rho)*np.maximum(np.abs(rho)-alpha,np.zeros((len(rho))))/z[ii]
            
            # 3)            
            # Add contribution to residual of updated coefficient
            # R -= x[ii]*A[:,ii]
            R -= np.dot(A[:,ii*num_maxima:(ii+1)*num_maxima],x[ii*num_maxima:(ii+1)*num_maxima])
#            plt.figure(1)            
#            plt.plot(R)
#            plt.pause(0.25)
            # update the maximum absolute coefficient update
            d_x_ii = np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima] - x_old))
            if d_x_ii > d_x_max:
                d_x_max = d_x_ii

            if np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima])) > x_max:
                x_max = np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima]))

#            x_pad = np.zeros(m)
#            x_pad[maxima] = x[ii*num_maxima:(ii+1)*num_maxima]
#            b_hat = np.zeros(m)
#            for i in range(num_vars):
#                b_hat += np.fft.ifft(np.multiply(Aft[:,i],np.fft.fft(x_pad))).real
#            
#            plt.figure(3)
#            plt.plot(b_hat)
#            plt.pause(0.25)
#        
        if (x_max == 0.0 or
            d_x_max / x_max < d_x_tol or
            iters == max_iter - 1):
            # the biggest coordinate update of this iteration was smaller
            # than the tolerance: check the duality gap as ultimate
            # stopping criterion

            # AtR = np.dot(A.T, R)
            AtR = np.dot(A_tran, R)


            if positive:
                dual_norm_AtR = np.max(AtR)
            else:
                dual_norm_AtR = np.max(np.abs(AtR))

            R_norm2 = np.dot(R, R)

            x_norm2 = np.dot(x, x)

            if (dual_norm_AtR > alpha):
                const = alpha / dual_norm_AtR
                A_norm2 = R_norm2 * (const ** 2)
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            l1_norm = np.sum(x)

            # np.dot(R.T, b)
            gap += (alpha * l1_norm - const * np.dot(R.T,b))

            if gap < tol:
                break
                # return if we reached desired tolerance
                
    return x

def lasso_solver(A,
                 A_tran,
                 A0,
                 A0_tran,
                 maxima,
                 b,
                 alpha, 
                 max_iter,
                 positive=0):

    m = A.shape[0]
    n = A.shape[1]
    num_vars = A0.shape[1]
    num_maxima = maxima.shape[0]

    # precompute sum of sqaures of columns
    z = np.sum(A**2,0) 
   
    # initialize x
    x = np.zeros(n)
    x_old = np.zeros(num_maxima)       
    x_pad = np.zeros(m)        
    AtR = np.zeros(m)
    R = np.zeros(m) 
    Rft = np.zeros(m) 


    rho = np.zeros(num_maxima)   
    tol = 1e-4
    gap = tol + 1.0
    d_x_tol = tol  
    

    Aft = np.fft.fft(A0,axis=0)
    A0_tran_ft = np.fft.fft(A0_tran,axis=0)       
       
    # R = b - np.dot(A,x) but x is zeros   
    R = b.copy()
    Rft = b.copy()
    
    # tol *= np.dot(b,b)
    tol *= np.dot(b,b)
 
    for iters in range(max_iter):
        x_max = 0.0
        d_x_max = 0.0
        
        for ii in range(num_vars):
            # 1)            
            x_old = x[ii*num_maxima:(ii+1)*num_maxima]
            x_pad = np.zeros(m)
            x_pad[maxima] = x_old
            x_ft = np.fft.fft(x_pad)
            
            # Remove contribution to residual of current coefficient
            # R += x_ii*A[:,ii]
            R += np.fft.ifft(np.multiply(Aft[:,ii],x_ft)).real
#            plt.figure(1)
#            plt.plot(R)
#            plt.pause(0.25)
            # 2)
            # Soft threshold update with optional positivity constraint
            rho = np.fft.ifft(np.multiply(A0_tran_ft[:,ii],np.fft.fft(R))).real[maxima]          
#            plt.figure(2)
#            plt.plot(rho,'o')
#            plt.pause(0.25)            
            # rho = np.dot(A[:,ii],R)
            if positive:
                rho[rho<0]=0
            x[num_maxima*ii:num_maxima*(ii+1)] = np.sign(rho)*np.maximum(np.abs(rho)-alpha,np.zeros((len(rho))))/z[ii]
            
            # 3)
            x_pad = np.zeros(m)
            x_pad[maxima] = x[ii*num_maxima:(ii+1)*num_maxima]
            x_ft = np.fft.fft(x_pad)
            
            # Add contribution to residual of updated coefficient
            # R -= x[ii]*A[:,ii]
            R -= np.fft.ifft(np.multiply(Aft[:,ii],x_ft)).real
#            plt.figure(1)            
#            plt.plot(R)
#            plt.pause(0.25)
            # update the maximum absolute coefficient update
            d_x_ii = np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima] - x_old))
            if d_x_ii > d_x_max:
                d_x_max = d_x_ii

            if np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima])) > x_max:
                x_max = np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima]))

#            x_pad = np.zeros(m)
#            x_pad[maxima] = x[ii*num_maxima:(ii+1)*num_maxima]
#            b_hat = np.zeros(m)
#            for i in range(num_vars):
#                b_hat += np.fft.ifft(np.multiply(Aft[:,i],np.fft.fft(x_pad))).real
#            
#            plt.figure(3)
#            plt.plot(b_hat)
#            plt.pause(0.25)
#        
        if (x_max == 0.0 or
            d_x_max / x_max < d_x_tol or
            iters == max_iter - 1):
            # the biggest coordinate update of this iteration was smaller
            # than the tolerance: check the duality gap as ultimate
            # stopping criterion

            # AtR = np.dot(A.T, R)
            AtR = np.zeros(m)
            Rft = np.fft.fft(R)
            for i in range(num_vars):
                AtR += np.real(np.fft.ifft(np.multiply(A0_tran_ft[:,i],Rft)))


            if positive:
                dual_norm_AtR = np.max(AtR)
            else:
                dual_norm_AtR = np.max(np.abs(AtR))

            R_norm2 = np.dot(R, R)

            x_norm2 = np.dot(x, x)

            if (dual_norm_AtR > alpha):
                const = alpha / dual_norm_AtR
                A_norm2 = R_norm2 * (const ** 2)
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            l1_norm = np.sum(x)

            # np.dot(R.T, b)
            gap += (alpha * l1_norm - const * np.dot(R.T,b))

            if gap < tol:
                break
                # return if we reached desired tolerance
                
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
max_iters = 200

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


print('py_time')
print(py_time)
#print('cy_time')
#print(cy_time)
print('sci_time')
print(sci_time)
print('batch_time')
print(batch_time)


#%%
#norm(B.dot(lasso.coef_) - signal)/norm(signal)
#norm(B.dot(coefs1) - signal)/norm(signal)
#norm(B.dot(coefs2) - signal)/norm(signal)

plt.plot(signal,'-og')
#plt.plot(y_hat1,'-ob')

plt.plot(y_hat3,'-or')
plt.plot(y_hat4,'-sk')
plt.plot(y_hat2,'-xm')

