# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 13:35:59 2017

@author: dbanco02
"""

## Init
from __future__ import division
from __future__ import print_function

import numpy as np
from numpy.linalg import norm 
import matplotlib.pyplot as plt
import os
import time

from scipy.signal import argrelmax
from sklearn.linear_model import Lasso

import RingImageProcessing as RingIP
import EllipticModels as EM
import DataAnalysis as DA
import LassoSolvers


#%% Create synthetic test data

# Create cartesian test data
num_r = 500
num_theta = 500

# Azimuthal units
theta = np.linspace(0,2*np.pi,num_theta)
dtheta = theta[1]-theta[0]

# Pixel units
r = np.linspace(0,num_r-1,num_r)
dr = 1

img = np.zeros((num_r,num_theta))

# Generate coordiantes for each pixel
r_coord = r.copy()
for i in range(num_theta-1):
    r_coord = np.vstack((r_coord,r))
r_coord = np.transpose(r_coord)

theta_coord = theta.copy()
for i in range(num_r-1):
    theta_coord = np.vstack((theta_coord,theta))
   
#%% Stack up a bunch of test images shiftes in theta using 2D gaussian basis
m_r = num_r/2
m_theta = np.linspace(0,num_theta,num_r)
v_r =     num_r/2
v_theta = num_theta/2
images = np.zeros((num_theta,0))
for j in range(num_r):
    dr = 1
    dtheta = 1
    B = EM.gaussian_basis_wrap_2D(num_theta,dtheta,m_theta[j], v_theta, 
                                  num_r,    dr,    m_r,        v_r)
    img = B.reshape((num_theta,num_r))
    images = np.append(images,img,1)    
plt.imshow(images)

#%% Test 2D convolution to generate 2D shifted Gaussian basis functions
num_r = 500
num_theta = 500
v_r =     num_r/2
v_theta = num_theta/2
m_r = 0
m_theta = 0
dr = 1
dtheta = 1
B0 = EM.gaussian_basis_wrap_2D(num_theta,dtheta,m_theta,    v_theta, 
                               num_r,    dr,    m_r,        v_r).reshape((num_theta,num_r))

delta2D = np.zeros((num_theta,num_r))
for i in np.arange(50,450,50):
    delta2D[i,i] = 1
         
result = LassoSolvers.Ax_ft_2D(B0,delta2D)
plt.imshow(np.abs(result))
# Success


#%% Test Ax and AtR on small test data
# Stack up a bunch of test images shiftes in theta using 2D gaussian basis

num_theta= 10
dtheta = 2*np.pi/num_theta

dr = 5
num_rad = 2*dr
drad = 1

# Create fitting matrix 
reload(EM)
reload(LassoSolvers)
num_var_t = 5
num_var_r = 5
var_theta = np.linspace((dtheta),(np.pi/16),num_var_t)**2
var_rad   = np.linspace(drad,10,num_var_r)**2
                       
#Store stack of image slices
B0_stack = np.zeros((num_rad,num_theta,num_var_t,num_var_r))
B0_stack_ft = np.zeros((num_rad,num_theta,num_var_t,num_var_r))
B_full = np.zeros((num_rad,num_theta,num_var_t,num_var_t,num_rad,num_theta))
m_r = 0
m_theta = 0
for t in range(num_var_t):
    for r in range(num_var_r):
        B0 = EM.gaussian_basis_wrap_2D(num_theta,dtheta,  m_theta,    var_theta[t], 
                                       num_rad,  drad,    m_r,        var_rad[r]).reshape((num_rad,num_theta))
        B0_stack[:,:,t,r] = B0
        B0_stack_ft[:,:,t,r] = np.fft.fft2(B0)
        for i in range(num_rad):
            for j in range(num_theta):
                B0_shift = np.concatenate([B0[-i:,:], B0[0:-i,:]],axis=0)
                B_shift = np.concatenate([B0_shift[:,-j:], B0_shift[:,0:-j]],axis=1)
                B_full[:,:,t,r,i,j] = B_shift    


B0ft_list = []
for tv in range(num_var_t):
    for rv in range(num_var_r):
        B0ft_list.append([np.fft.fft2(B0_stack[:,:,tv,rv]).real])

x = np.random.rand(num_rad,num_theta,num_var_t,num_var_t)
R = np.random.rand(num_rad,num_theta)
Ax = np.zeros((num_rad,num_theta))
AtR = np.zeros((num_rad,num_theta,num_var_t,num_var_t))

# Compute Ax and AtR regularly
for t in range(num_var_t):
    for r in range(num_var_r):
        for i in range(num_rad):
            for j in range(num_theta):      
                Ax += B_full[:,:,t,r,i,j]*x[i,j,t,r]
                AtR[i,j,t,r] = np.sum(B_full[:,:,t,r,i,j]*R)
                
# Compute Ax and AtR with convolution
Ax_conv = LassoSolvers.Ax_ft_2D_Parallel(B0ft_list,x)
AtR_conv = LassoSolvers.AtR_ft_2D_Parallel(B0ft_list,R)

print('Ax Error: '  + str(norm(Ax-Ax_conv)/norm(Ax)))
print('AtR Error: ' + str(norm(AtR-AtR_conv)/norm(AtR)))
#%% Attempt to obtain first pixel of every shifted basis
reload(EM)
B0_regular = EM.gaussian_basis_wrap_2D(num_theta,dtheta,  m_theta,    var_theta[14], 
                                       num_rad,  drad,    m_r,        var_rad[9]).reshape((num_rad,num_theta))
B0_tran_tube, domain = EM.gaussian_basis_wrap_2D_tube(num_theta,
                                                    dtheta,
                                                    var_theta[14],
                                                    num_rad,
                                                    drad,
                                                    var_rad[9])

B0_tran_slice = B0_tran_tube.reshape((num_rad,num_theta))
domain_theta = domain[0].reshape((num_rad,num_theta))
domain_rad = domain[1].reshape((num_rad,num_theta))

print(norm(B0_regular-B0_tran_slice))
print(np.argmax(B0_regular-B0_tran_slice))

plt.figure(1)
plt.imshow(B0_tran_slice,cmap='jet', interpolation='nearest')

plt.figure(2)
plt.imshow(domain_theta,cmap='jet', interpolation='nearest')
    
plt.figure(3)
plt.imshow(domain_rad,cmap='jet', interpolation='nearest')      

plt.figure(4)
plt.imshow(B0_regular,cmap='jet', interpolation='nearest')

plt.figure(5)
plt.imshow(B0_regular-B0_tran_slice,cmap='jet', interpolation='nearest')

# Success

# Test to make sure B0 and B0_tran are identical
for t in range(num_theta):
    for r in range(num_rad):
        B0_regular = EM.gaussian_basis_wrap_2D(num_theta,dtheta,  t,    var_theta[14], 
                                               num_rad,  drad,    r,        var_rad[9])
        B0_tran_tube, domain = EM.gaussian_basis_wrap_2D_shift_tube(num_theta,
                                                            dtheta,
                                                            var_theta[14],
                                                            num_rad,
                                                            drad,
                                                            var_rad[9],
                                                            t, 
                                                            r)
        print(norm(B0_regular-B0_tran_tube))
        
# Succes


#%% Test low rank decomposition of gaussian basis functions
# Create fitting matrix 
reload(EM)
num_var_t = 15
num_var_r = 10
var_theta = np.linspace((dtheta),(np.pi/32),num_var_t)**2
var_rad   = np.linspace(drad,3,num_var_r)**2

# Generate unshifted basis function for each variance combination

#Store stack of image slices
B0_stack = np.zeros((num_rad,num_theta,num_var_t,num_var_r))
for t in range(num_var_t):
    for r in range(num_var_r):
        m_r = 0
        m_theta = 0
        B0 = EM.gaussian_basis_wrap_2D(num_theta,dtheta,  m_theta,    var_theta[t], 
                                       num_rad,  drad,    m_r,        var_rad[r]).reshape((num_rad,num_theta))
        B0_stack[:,:,t,r] = B0       

B0_stack_ft = np.zeros(B0_stack.shape)
for tv in range(num_var_t):
    for rv in range(num_var_r):
        B0_stack_ft[:,:,tv,rv] = np.fft.fft2(B0_stack[:,:,tv,rv]) 

u,s,v = np.linalg.svd(np.fft.fft2(B0_stack[:,:,14,9]))

u[:,1:] = 0
s[1:] = 0
v[1:,:] = 0
 
S = np.zeros((u.shape[0],v.shape[0]))
 
norm(B0_stack[:,:,14,9] - np.dot(u,S).dot(v.T) )/norm(B0_stack[:,:,14,9])

plt.figure(3)
plt.plot(s,'o')

# Conclusion: basis function is rank 1
# Doesn't seem like there are computational savings, but there are storage
# savings that could facilitate parallel processing