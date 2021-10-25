# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:56:01 2016

@author: Dan
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

### Data directory
data_dir = os.path.join('E:','CHESS_raw_data')
out_dir  = os.path.join(data_dir,'out')

#%% Load Real Image Data

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
radius                = 718                # ring radius in pixels
dr                    = 30                 # half of ring width in pixels
min_amp               = 25                 # minimum acceptable peak amplitude
vec_frac              = 0.25               # fraction of peaks that must be acceptable

sample    = DA.Specimen(specimen_name, data_dir, out_dir,
                        step_names, dic_files, dark_dirs, 
                        init_dirs, dic_center, 
                        x_range,x_num, y_range, y_num,
                        detector_dist, true_center, 
                        e_rng, p_rng, t_rng, E, G, v)   

# Data path
struc_data = os.path.join('E:','CHESS_data','results_exp_aluminum_fixed.npy')
rsf_out = np.load(struc_data)

# Specify image data
img_num = 165
load_i = 3
rsf_ex = rsf_out[load_i][img_num]

# 1D azimuthal diffraction data
b = rsf_ex.f

# Parameters
num_theta = len(b)
dtheta = 2*np.pi/num_theta


#%% Generate Data instead
"""
num_theta = 300
dtheta = 2*np.pi/num_theta

scale = [200, 400, 330]
means = [40, 140, 200]
stds = [5, 10, 7]

np.random.normal(0,1,num_theta)

noise = np.random.normal(0,1,num_theta)
b = np.zeros(num_theta)

for i in range(3):
    b += scale[i]*EM.gaussian_basis_wrap(num_theta,
                                           dtheta,
                                           means[i],
                                           (stds[i]*dtheta)**2)
"""
#%% Lasso fitting for each ring
num_var = 15

# Extract local maxima
maxima =argrelmax(b)[0]

# Create Gaussian basis function fitting matrix
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

B = np.array(B, order='F', copy=True)
B_tran = np.array(B_tran, order='F', copy=True)
B0 = np.array(B0, order='F',copy=True)
B0_tran = np.array(B0_tran, order='F',copy=True)

# Fourier transform of matrices
#B0ft = np.fft.fft(B0,axis=0)
#B0_tranft = np.fft.fft(B0_tran,axis=0)
#B0ft = np.array(B0ft, order='F',copy=True)
#B0_tranft = np.array(B0_tranft, order='F',copy=True)   

# Output as .mat files for MATLAB
#from scipy.io import savemat
#savemat('test_mats_small.mat',dict(B=B, b=b))



#%% Test Lasso fitting routines
reload(RingIP)
reload(LassoSolvers)

# Test outputs
num_tests = 5
times = np.zeros((num_tests))
rel_error = np.zeros((num_tests))
coefs = np.zeros((B.shape[1],num_tests))
y = np.zeros((num_theta,num_tests))

# Parameters
l1_ratio = 0.08
max_iters = 800

test_names = ['Coordinate Ascent',
              'FISTA',
              'Circulant FISTA',
              'Scikit Coordinate Ascent',
              'SciPy Minimize']

reload(LassoSolvers)
print('Coordinate Ascent')
coefs1, times1 = LassoSolvers.cython_coord_ascent(B,
                                              b,
                                              l1_ratio, 
                                              max_iters,
                                              positive=1, 
                                              benchmark=1)
coefs[:,0] = coefs1
times[0]   = times1


print('FISTA')
coefs1, times1 = LassoSolvers.fista(B, b,
                             l1_ratio, max_iters,  
                             positive=1, benchmark=1)
coefs[:,1] = coefs1
times[1]   = times1


print('Circulant FISTA')
coefs1, times1 = LassoSolvers.fista_circulant(B, B0, B0_tran, maxima, b, 
                                       l1_ratio, max_iters, 
                                       positive=1, benchmark=1) 
coefs[:,2] = coefs1
times[2]   = times1

print('Scikit Coordinate Ascent')
t1 = time.time()
lasso = Lasso(alpha=l1_ratio,
              max_iter=max_iters,
              fit_intercept=0,
              positive=True)
lasso.fit(B,b)
y_hat3 = B.dot(lasso.coef_)
t2 = time.time()  
sci_time = t2-t1
coefs[:,3] = lasso.coef_
times[3]   = sci_time

reload(RingIP)
reload(LassoSolvers)
print('SciPy Minimize')

eps = 10**(-8)
coefs1, times1 = LassoSolvers.LASSO_approx_tnc(B0,B0_tran,maxima,
                                               num_var,num_theta,
                                               b,l1_ratio,eps,
                                               max_iters,
                                               callback=1,
                                               benchmark=1)
coefs[:,4] = coefs1[0]
times[4]   = times1
obj_hist = np.load('obj_hist.npy')
plt.plot(obj_hist)
     # plot cost as func of iter
     

reload(RingIP)
reload(LassoSolvers)

for i in range(num_tests):
    y[:,i] = np.dot(B,coefs[:,i])
    rel_error[i] = norm(y[:,i]-b)/norm(b)
    
#%% Plot Stuff

for i in range(num_tests):
    plt.figure(i)
    plt.plot(b,'-ob')
    plt.plot(y[:,i],'-or')
    plt.title(test_names[i])
    
    
    
#%%
"""
import scipy.optimize as opt
x0 = np.zeros(num_var*len(maxima))
m = num_theta
n = B0.shape[1]
A0_ft = np.fft.fft(B0,axis=0)    
A0_tran_ft = np.fft.fft(B0_tran,axis=0)  

# Construct positivity bounds
bnd = []
for i in range(len(x0)):
    bnd.append((0,None))

x = opt.fmin_l_bfgs_b(LassoSolvers.obj_func_circ,
                       x0,
                       args=(A0_ft,A0_tran_ft,b,maxima,num_var,m,n,l1_ratio,eps),
                       bounds=bnd,
                       maxiter=max_iters,
                       disp=1,
                       iprint=1)

 """
