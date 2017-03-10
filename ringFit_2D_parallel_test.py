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


#%% Test Ax and AtR on small test data
# Stack up a bunch of test images shiftes in theta using 2D gaussian basis

num_theta= 500
dtheta = 2*np.pi/num_theta

dr = 5
num_rad = 2*dr
drad = 1

# Create fitting matrix 
reload(EM)
reload(LassoSolvers)
num_var_t = 15
num_var_r = 10
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

x = np.random.rand(num_rad,num_theta,num_var_t,num_var_r)
R = np.random.rand(num_rad,num_theta)
Ax = np.zeros((num_rad,num_theta))
AtR = np.zeros((num_rad,num_theta,num_var_t,num_var_r))
AtR_conv_para = np.zeros((num_rad,num_theta,num_var_t,num_var_r))

# Compute Ax and AtR regularly
for t in range(num_var_t):
    for r in range(num_var_r):
        for i in range(num_rad):
            for j in range(num_theta):      
                Ax += B_full[:,:,t,r,i,j]*x[i,j,t,r]
                AtR[i,j,t,r] = np.sum(B_full[:,:,t,r,i,j]*R)
                
# Compute Ax and AtR with convolution
a = time.time()
AtR_conv = LassoSolvers.AtR_ft_2D(B0_stack_ft,R)
timeAtR = time.time() - a

a = time.time()
Ax_conv = LassoSolvers.Ax_ft_2D(B0_stack_ft,x)
timeAx = time.time() - a

a = time.time()
AtR_conv_para2 = LassoSolvers.AtR_ft_2D_Parallel(B0ft_list[:],R)
timeAtR_para = time.time() - a

a = time.time()
Ax_conv_para = LassoSolvers.Ax_ft_2D_Parallel(B0ft_list[:],x)
timeAx_para = time.time() - a


k = 0
for tv in range(num_var_t):
    for rv in range(num_var_r):
    
        AtR_conv_para[:,:,tv,rv] = AtR_conv_para2[k,0,:,:]
        k += 1
    
print('Ax Error: '  + str(norm(Ax-Ax_conv)/norm(Ax)))
print('AtR Error: ' + str(norm(AtR-AtR_conv)/norm(AtR)))
print('Ax Parallel Error: '  + str(norm(Ax-Ax_conv_para)/norm(Ax)))
print('AtR Parallel Error: ' + str(norm(AtR-AtR_conv_para)/norm(AtR)))

print('Ax Time: '  + str(timeAx))
print('AtR Time: ' + str(timeAtR))
print('Ax Parallel Time: '  + str(timeAx_para))
print('AtR Parallel Time: ' + str(timeAtR_para))