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
import multiprocessing
from functools import partial

# Data, Interpolation, Fitting Parameters
load_step = 4
img_nums = range(0,50)

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

# Compute basis matrix and interpolate ring to polar coordinates 
ringModel = RM.RingModel(load_step, img_nums[0], radius-dr, radius+dr, 
		              num_theta, num_rad, dtheta, drad, var_theta, var_rad)
A0ft_stack = EM.unshifted_basis_matrix_ft_stack(ringModel.var_theta,
						                               ringModel.var_rad,
						                               ringModel.dtheta,
						                               ringModel.drad,
						                               ringModel.num_theta, ringModel.num_rad)

ringModel.l1_ratio = 1
ringModel.max_iters = 500
ringModel.lipschitz = 6e5

out_path = os.path.join('..','CHESS_results','ringModel_out')

ringModel_list = []
for img_num in img_nums:
    img_file = 'polar_image_al7075_load_'+str(load_step)+'_img_'+str(img_num)+'.npy'
    img_path = os.path.join('..','CHESS_data',img_file)
    ringModel = RM.RingModel(load_step, img_num, radius-dr, radius+dr, 
    		              num_theta, num_rad, dtheta, drad, var_theta, var_rad)
    ringModel.polar_image = np.load(img_path)
    ringModel.l1_ratio = 1
    ringModel.max_iters = 500
    ringModel.lipschitz = 6e5
    ringModel_list.append(ringModel)

partial_fit = partial(RM.fit_circulant_FISTA_Multiprocess,A0ft_stack=A0ft_stack,out_path=out_path,
positive=1,benchmark=1,verbose=2)

pool = multiprocessing.Pool()
pool.map(partial_fit,ringModel_list)

