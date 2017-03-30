# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 13:35:59 2017

@author: dbanco02
"""

## Init
from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplt as plt
import os
import RingModel as RM
import EllipticModels as EM
import multiprocessing
from functools import partial

# Data, Interpolation, Fitting Parameters
dr = 30
radius = 370
 
num_theta= 2048
num_rad = 2*dr

num_var_t = 15
num_var_r = 10

dtheta = 2*np.pi/num_theta
drad = 1

var_theta = np.linspace((dtheta),(np.pi/32),num_var_t)**2
var_rad   = np.linspace(drad,3,num_var_r)**2

A0_stack = EM.unshifted_basis_matrix_stack(var_theta,
                                           var_rad,
                                           dtheta,
                                           drad,
                                           num_theta, 
                                           num_rad)

result_path = os.path.join('..','CHESS_results')

var_signal = np.zeros((num_var_r,num_var_t,5,41,5))
for i, step in enumerate(range(5)):
    for j, img_num in enumerate(range(205)):
        file_name = 'ringModel_load_' + str(step)+'_img_' + str(img_num) + '.npy'
        file_path = os.path.join(result_path,file_name)
        ringModel = np.load(file_path)
        var_signal[:,:,i,j//5,j%5] = np.sum(np.sum(ringModel.coefs,0),0)*np.sum(np.sum(A0_stack,0),0)
        
        
high_var_rad = np.sum(np.sum(var_signal[3::,:,:,:,:],0),0)
high_var_theta = np.sum(np.sum(var_signal[:,3::,:,:,:],0),0)
mu = np.mean(high_var_rad[:,2:41,:].ravel())
sig = np.std(high_var_rad[:,2:41,:].ravel())
plt.figure(3)
plt.subplot(1,5, 1)  
plt.imshow(high_var_rad[i], vmin=0,vmax=1, interpolation='nearest')
plt.title('HighVar')
plt.axis('off')
plt.colorbar()    