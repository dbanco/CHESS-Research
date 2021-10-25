# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:01:01 2016

Batch peak fitting

@author: Dan
"""

## Init
from __future__ import division

import numpy as np
from numpy.linalg import norm 
import matplotlib.pyplot as plt
import os
import time

from scipy.signal import argrelmax
from sklearn.linear_model import Lasso

import DataReader
import RingImageProcessing as RingIP
import EllipticModels as EM
import PeakFitting as peak
import DataAnalysis as DA
reload(RingIP)
reload(EM)

### Data directory
data_dir = os.path.join('F:','CHESS')
out_dir  = os.path.join(data_dir,'out')

#%%
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

#%% Lasso fitting for each ring
# Specify image data
load_step_list = [3]
img_num_list = range(23,24)

num_theta = 2400
dtheta = 2*np.pi/num_theta
num_var = 100

var_domain = np.linspace((dtheta),(np.pi/16),num_var)**2

basis_path = os.path.join('basis_radial_sum','gaus_basis_shift_')

# Load one-dimensional fits
rsf_out = np.load('results_exp_final_radial_sum_fixed.npy')

for load_i in load_step_list:
    for img_num in img_num_list:
        print('Load step ' + str(load_i) +', img num ' + str(img_num) )
        
        rsf = rsf_out[load_i][img_num]        
        # Load ring image
        ring_image = sample.load_image(img_num,load_i) 
        
        # Determine support of ring from 1d ellipse fit
        distance = 8
        radius_low = min(rsf.a,rsf.b) - distance
        radius_high = max(rsf.a,rsf.b) + distance
        mask_image = RingIP.distance_threshold(ring_image,
                                               radius_low,
                                               radius_high,
                                               center=[rsf.xc,rsf.zc])
        x_mask,y_mask,f_mask = RingIP.get_points(ring_image,mask_image)
        num_samples = len(f_mask)        
        
        # Extract means and variances of fitted basis functions
        keep = rsf.coefs > 0
        az_means = rsf.means[keep]
        az_variances = rsf.variances[keep]
        az_coefs = rsf.coefs[keep]
        az_radii = np.sqrt((rsf.a*np.cos(az_means))**2 + 
                           (rsf.b*np.sin(az_means))**2)
        az_x0 = rsf.a*np.cos(az_means) + rsf.xc
        az_y0 = rsf.b*np.sin(az_means) + rsf.zc
        
        # Define domain of radial means and variances 
        num_rad_var = 5
        rad_variance_domain = np.linspace(1,7,num_rad_var)**2 
        rad_distances = np.array([1,2,3])
        B = np.empty([num_samples,0])
        B_az_mean = np.empty(0)
        B_az_var  = np.empty(0)
        B_rad_dx = np.empty(0)
        B_rad_dy = np.empty(0)
        B_rad_var  = np.empty(0)
        
        for i_az, az in enumerate(az_means):       
            print(i_az)
            x_radii, y_radii = RingIP.line_normal_to_ellipse(az_x0[i_az],
                                                             az_y0[i_az],
                                                             rsf.a,rsf.b,
                                                             rad_distances)
            dx = az_x0[i_az] - x_radii.ravel()
            dy = az_y0[i_az] - y_radii.ravel()
            for j in range(len(dx)):
                for rad_var in rad_variance_domain:
                    B0 = EM.ellipse_basis_single_spot([x_mask,y_mask],rsf.a,rsf.b,
                                                 rsf.xc,rsf.zc,dx[j],dy[j],
                                                 az,num_theta,rad_var,
                                                 az_variances[i_az])
                    B = np.hstack([B,B0[:,None]])
                    B_az_mean = np.append(B_az_mean,az)
                    B_az_var  = np.append(B_az_var,az_variances[i_az])
                    B_rad_dx = np.append(B_rad_dx,dx[j])
                    B_rad_dy = np.append(B_rad_dy,dy[j])
                    B_rad_var  = np.append(B_rad_var,rad_var)

        #wsqanp.save(os.path.join(name,'full_ellipse_basis.npy'),B)
        # might need to adjust the units of az (azimuthal mean)                                         
                                         
        # Construct 2d elliptic basis matrix based on 1d fit and define the
        # basis path. Only `
        
        # Try and make it stored as a sparse matrix? How can I tell if this
        # will be advantageous?????
        r2Df = RingIP.Ring2DFit(load_i,img_num,rsf.a,rsf.b,rsf.xc,rsf.zc,
                    f_mask,[x_mask, y_mask],ring_image.shape,
                    B_az_mean,B_az_var,B_rad_dx,
                    B_rad_dy,B_rad_var,0,0,0,0,
                    0,0,0)
               
        r2Df2,_ = RingIP.circulant_lasso_ring_2Dfit(r2Df,B)
           
        #r2Df_out[load_i].append(r2Df2)

    #np.save('initial_test_elliptic_2Dfit.npy',r2Df2)

# Messin' around
reload(RingIP)
reload(EM)
r2Df_test = RingIP.Ring2DFit(r2Df2.load_step,
                        r2Df2.img_num,
                        r2Df2.a,
                        r2Df2.b,
                        r2Df2.xc,
                        r2Df2.yc,
                        r2Df2.f,
                        r2Df2.mask_xy,
                        r2Df2.img_shape,
                        r2Df2.az_means,
                        r2Df2.az_variances,
                        r2Df2.rad_dx,
                        r2Df2.rad_dy,
                        r2Df2.rad_variances,
                        r2Df2.coefs,
                        r2Df2.fit,
                        r2Df2.l1_ratio,
                        r2Df2.n_iter,
                        r2Df2.fit_error,
                        r2Df2.rel_fit_error)
r2Df3.display_fit_image(0)
r2Df3.display_f_image(1)

az_var_domain = np.linspace((dtheta),(np.pi/16),num_var)**2
rad_var_domain = np.linspace(1,7,num_rad_var)**2 

#%%
reload(RingIP)
reload(EM)
B0 = 400*EM.ellipse_basis_single_spot(r2Df_test.mask_xy,r2Df_test.a,r2Df_test.b,
                                                 r2Df_test.xc,r2Df_test.yc,
                                                 0,0,
                                                 0,2400,
                                                 3,
                                                 az_var_domain[90])     
r2Df_test.fit = B0           
r2Df_test.display_fit_image(0)
plt.xlim([np.min(950),np.max(1450)])
plt.ylim([np.min(950),np.max(1450)])
r2Df_test.display_f_image(1)
plt.xlim([np.min(950),np.max(1450)])
plt.ylim([np.min(950),np.max(1450)])

#dx[j],
#dy[j],
#az,num_theta,
#rad_var,
#az_variances[i_az]
                                                 
