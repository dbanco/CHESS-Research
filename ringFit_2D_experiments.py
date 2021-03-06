# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 13:35:59 2017

@author: dbanco02
"""

## Init
from __future__ import division
from __future__ import print_function

import numpy as np
import os
import RingModel as RM
import EllipticModels as EM



# Data, Interpolation, Fitting Parameters
load_step = 0
img_num = 35

img_file = 'polar_image_al7075_load_'+str(load_step)+'_img_'+str(img_num)+'.npy'
img_path = os.path.join('/home','dbanco','CHESS_data',img_file)
 
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

# Initialize Ring Model
ringModel = RM.RingModel(load_step, img_num, radius-dr, radius+dr, 
    		              num_theta, num_rad, dtheta, drad, var_theta, var_rad)

# Compute basis matrix and interpolate ring to polar coordinates 
A0_stack = EM.unshifted_basis_matrix_stack(ringModel.var_theta,
						                               ringModel.var_rad,
						                               ringModel.dtheta,
						                               ringModel.drad,
						                               ringModel.num_theta, ringModel.num_rad)

A0ft_svd_list = ringModel.generate_basis_matrices()

ringModel.polar_image = np.load(img_path)
print(ringModel.polar_image.shape)
ringModel.l1_ratio = 1
ringModel.max_iters = 500
ringModel.compute_lipschitz(A0_stack)

ringModel.fit_circulant_FISTA(A0ft_svd_list,positive=1,benchmark=1,verbose=1)

np.save('ringModel_0_35.npy',ringModel)

