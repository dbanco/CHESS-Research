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
import CirculantOperations as CO

#%% Test Ax and AtR on small test data
# Stack up a bunch of test images shiftes in theta using 2D gaussian basis

num_theta= 25
dtheta = 2*np.pi/num_theta

dr = 6
num_rad = 2*dr
drad = 1

# Create fitting matrix 
reload(EM)
reload(LassoSolvers)
num_var_t = 8
num_var_r = 6
var_theta = np.linspace((dtheta),(np.pi/16),num_var_t)**2
var_rad   = np.linspace(drad,10,num_var_r)**2
                       
#Store stack of image slices
B0ft_stack    = EM.unshifted_basis_matrix_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad) 
B0ft_list     = EM.unshifted_basis_matrix_list(var_theta,var_rad,dtheta,drad,num_theta,num_rad)
B0ft_list_svd = EM.unshifted_basis_svd_list(var_theta,var_rad,dtheta,drad,num_theta,num_rad)


x = np.random.rand(num_rad,num_theta,num_var_t,num_var_r)
R = np.random.rand(num_rad,num_theta)
Ax = np.zeros((num_rad,num_theta))
AtR = np.zeros((num_rad,num_theta,num_var_t,num_var_r))
AtR_conv_para = np.zeros((num_rad,num_theta,num_var_t,num_var_r))
     
# Convolution
a = time.time()
AtR_conv = CO.AtR_ft_2D(B0ft_stack,R)
timeAtR = time.time() - a

a = time.time()
Ax_conv = CO.Ax_ft_2D(B0ft_stack,x)
timeAx = time.time() - a

# Convolution + Parallel
a = time.time()
AtR_conv_para2 = CO.AtR_ft_2D_Parallel(B0ft_list[:],R)
timeAtR_para = time.time() - a

a = time.time()
Ax_conv_para = CO.Ax_ft_2D_Parallel(B0ft_list[:],x)
timeAx_para = time.time() - a

# Convolution + Parallel + SVD
a = time.time()
AtR_conv_para2 = CO.AtR_ft_2D_Parallel_SVD(B0ft_list_svd[:],R)
timeAtR_para = time.time() - a

a = time.time()
Ax_conv_para = CO.Ax_ft_2D_Parallel_SVD(B0ft_list_svd[:],x)
timeAx_para = time.time() - a

print(AtR_conv_para2.shape)
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
