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


reload(EM)
reload(LassoSolvers)

# Data directory
data_dir = os.path.join('E:','CHESS_raw_data')
out_dir  = os.path.join(data_dir,'out')

# Load Real Image Data
specimen_name         = 'al7075_mlf'
step_names            = ['initial',  '1turn',    '2turn',    '3turn',    'unload']
dic_files             = ['dic_4536', 'dic_4537', 'dic_4538', 'dic_4539', 'dic_4540']
dark_dirs             = [ 68,         277,        485,        692,        899]
init_dirs             = [ 68,         277,        485,        692,        899]
dic_center            = [0.16, 2.9]       # sample center (mm) in vic2d coordinates 
x_range, x_num        = [-6, 6], 5         # x-ray measurements, range in mm
y_range, y_num        = [-5, 5], 41        # x-ray measurements, range in mm
detector_dist         = 3289.95            # pixels
true_center           = [1020.67, 1024.61] # [row, column] of detector image center in pixels (shifted by 1 for python index)
e_rng                 = [-0.012, 0.012]    # elastic strain range
p_rng                 = [-0.024, 0.024]    # plastic strain range
t_rng                 = [-0.036, 0.036]    # total strain range
E, G, v               =  71.7, 26.9, 0.33  # elastic modulus (GPa), shear modulus (GPa), poisson's ratio

ring_name             = 'al_311'
radius                = 370                # ring radius in pixels
dr                    = 30                 # half of ring width in pixels
min_amp               = 25                 # minimum acceptable peak amplitude
vec_frac              = 0.25               # fraction of peaks that must be acceptable

sample    = DA.Specimen(specimen_name, data_dir, out_dir,
                        step_names, dic_files, dark_dirs, 
                        init_dirs, dic_center, 
                        x_range,x_num, y_range, y_num,
                        detector_dist, true_center, 
                        e_rng, p_rng, t_rng, E, G, v) 

num_theta= 2048
dtheta = 2*np.pi/num_theta

num_rad = 2*dr
drad = 1

img_num = 35
load_step = 0
image = sample.load_image(img_num,load_step)

polar_image = np.zeros((num_rad,num_theta))

for i, r in enumerate(np.arange(radius-dr,radius+dr,1)):
    fi, theta_domain = RingIP.azimuthal_projection(image,true_center,r,0,2*np.pi-dtheta,num_theta)
    polar_image[i,:] = fi

#plt.figure(1)
#plt.imshow(polar_image,cmap='jet',vmin=0,vmax=200, interpolation='nearest')
#plt.figure(2)
#plt.imshow(image,cmap='jet',vmin=0,vmax=200, interpolation='nearest')

# Success


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

# Fit 2D data
l1_ratio = 1
max_iters = 500

#plt.figure(1)
#plt.imshow(8000*B0_stack[:,:,12,9],cmap='jet', interpolation='nearest',vmin=0,vmax=200)

#summ = 0
#for t in range(num_var_t):
#    for r in range(num_var_r):
#        eigt = np.linalg.eig( np.dot(B0_stack[:,:,t,r].T,
#                                    B0_stack[:,:,t,r] ))
#        summ += np.max(eigt[0].real)
#        print(summ)
# Sum was 7782.123
        
eig = np.linalg.eig( np.dot(B0_stack.reshape((num_rad*num_var_t*num_var_r,num_theta)).T,
                            B0_stack.reshape((num_rad*num_var_t*num_var_r,num_theta)) ))
L = np.max(eig[0].real)*num_rad*num_theta/1200
          
print('Circulant FISTA 2D')
x_hat, times = LassoSolvers.fista_circulant_2D(B0_stack, polar_image, 
                                       L, l1_ratio, max_iters, 
                                       positive=1, benchmark=1,verbose=1) 

y_hat = LassoSolvers.Ax_ft_2D(B0_stack_ft,x_hat)

# Result Error
error = norm(y_hat-polar_image)/norm(polar_image)
print(error)

# Result Sparsity 
pos_coef = np.sum(x_hat>0)
tot_coef = len(x_hat.ravel())
sparsity0 = pos_coef/tot_coef
one_coef = np.sum(x_hat>1)
sparsity1 = one_coef/tot_coef
print('x > 0         |  x > 1             | total x')
print(str(pos_coef) + '       |  ' + str(one_coef) + '             | ' + str(tot_coef))
print(str(sparsity0) + ' |  ' + str(sparsity1))

plt.figure(1)
plt.imshow(y_hat,cmap='jet', interpolation='nearest',vmin=0,vmax=200)

plt.figure(2)
plt.imshow(polar_image,cmap='jet', interpolation='nearest',vmin=0,vmax=200)

data_out_dir = os.path.join('E:','CHESS_data')
np.save(os.path.join(data_out_dir,'result_2D_500_12.npy'),
        (x_hat, times, y_hat, polar_image,
         var_theta,var_rad,L,l1_ratio,max_iters))


#%% Load and view results of previously completed trials 
data_out_dir = os.path.join('E:','CHESS_data')
print_file = os.path.join(data_out_dir,'result_2D_800_7.npy')
x_hat, times, y_hat, polar_image, var_theta,var_rad,L,l1_ratio,max_iters = np.load(print_file)

y_hat  = LassoSolvers.Ax_ft_2D(B0_stack_ft,x_hat)

# Result Error
error = norm(y_hat-polar_image)/norm(polar_image)
print('Error: ' + str(error))

# Result Sparsity 
pos_coef = np.sum(x_hat>0)
tot_coef = len(x_hat.ravel())
sparsity0 = pos_coef/tot_coef
one_coef = np.sum(x_hat>1)
sparsity1 = one_coef/tot_coef
print('x > 0         |  x > 1             | total x')
print(str(pos_coef) + '       |  ' + str(one_coef) + '             | ' + str(tot_coef))
print(str(sparsity0) + ' |  ' + str(sparsity1))
print('L: ' +str(L))
print('l1_ratio: '+ str(l1_ratio))
print('max_iters: ' +str(max_iters))
print('var_theta: '+ str(var_theta[-1]))
print('var_rad: '+ str(var_rad[-1]))
print('time: '+ str(times))

plt.figure(1)
plt.imshow(y_hat,cmap='jet', interpolation='nearest',vmin=0,vmax=200)

plt.figure(2)
plt.imshow(polar_image,cmap='jet', interpolation='nearest',vmin=0,vmax=200)
