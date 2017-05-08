# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 13:35:59 2017

@author: dbanco02
"""

## Init
from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import os
import RingModel as RM
import EllipticModels as EM

#%% Data, Interpolation, Fitting Parameters
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

A0_sum = np.sum(np.sum(A0_stack,0),0)

result_path = os.path.join('E:','CHESS_results','full_2D')

#%% Load results from files
#var_signal = np.zeros((num_var_t,num_var_r,5,41,5))
#rel_error = np.zeros((5,41,5))

load_steps = 5
for i, step in enumerate(range(3,4)):
    for j, img_num in enumerate(range(156,157)):
        print('Load: ' + str(step) + '   Image: ' + str(img_num))
        file_name = 'ringModel_out_load_' + str(step)+'_img_' + str(img_num) + '.npy'
        file_path = os.path.join(result_path,file_name)
        ringModel = np.load(file_path)
        rel_error[step,img_num//5,img_num%5] = ringModel[()].rel_fit_error
        var_signal[:,:,step,img_num//5,img_num%5] = np.sum(np.sum(ringModel[()].coefs,0),0)*A0_sum

        
#%% Plot results

cutoff_t = 6
cutoff_r = 4
total_var = np.sum(np.sum(var_signal[:,:,0:load_steps,:,:],0),0)
for i in range(load_steps):
    plt.figure(1)
    high_var_theta = np.sum(np.sum(var_signal[cutoff_t::,:,0:load_steps,:,:],0),0)/total_var
    mu = np.mean(high_var_theta[0:load_steps,2:41,:].ravel())
    sig = np.std(high_var_theta[0:load_steps,2:41,:].ravel())
    plt.subplot(1,5, i+1)  
    plt.imshow(high_var_theta[i], vmin=0,vmax=np.max(high_var_theta[0:load_steps,2:41,:].ravel())/10, interpolation='nearest')
    if(i ==2 ):
        plt.title('Theta Spread')
    plt.axis('off')
    #if(i==4):
     #   plt.colorbar()    
    
    plt.figure(2)
    high_var_rad = np.sum(np.sum(var_signal[:,cutoff_r::,0:load_steps,:,:],0),0)/total_var
    mu = np.mean(high_var_rad[0:load_steps,2:41,:].ravel())
    sig = np.std(high_var_rad[0:load_steps,2:41,:].ravel())
    plt.subplot(1,5, i+1)  
    plt.imshow(high_var_rad[i], vmin=0,vmax=np.max(high_var_rad[0:load_steps,2:41,:].ravel()), interpolation='nearest')
    if(i ==2 ):
        plt.title('Radial Spread')
    plt.axis('off')
    #if(i==4):
     #   plt.colorbar()    
    
    plt.figure(3)
    mu = np.mean(rel_error[0:load_steps,2:41,:].ravel())
    sig = np.std(rel_error[0:load_steps,2:41,:].ravel())
    plt.subplot(1,5, i+1)  
    plt.imshow(rel_error[i], vmin=0,vmax=1, interpolation='nearest')
    if(i ==2 ):
        plt.title('Fit Error')
    plt.axis('off')
    #if(i==4):
     #   plt.colorbar()
     
#%% Load polar and fit image 

img_num//5,img_num%5

load_step = 5
#row = 4
#col = 2
#img_num = row*
print('Load: ' + str(step) + '   Image: ' + str(img_num))
file_name = 'ringModel_out_load_' + str(step)+'_img_' + str(img_num) + '.npy'
file_path = os.path.join(result_path,file_name)
ringModel = np.load(file_path)

plt.figure(4)
plt.subplot(3,1,1) 
plt.imshow(ringModel[()].polar_image, vmin=0,vmax=200, interpolation='nearest')
plt.subplot(3,1,2)
plt.title(str(ringModel[()].rel_fit_error))
plt.imshow(ringModel[()].fit_image,   vmin=0,vmax=200, interpolation='nearest')

plt.subplot(3,1,3)
plt.imshow(500*A0_stack[:,:,14,9],   vmin=0,vmax=200, interpolation='nearest')
