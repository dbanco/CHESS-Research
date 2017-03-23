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
import random
from functools import partial
import multiprocessing

# Data, Interpolation, Fitting Parameters
load_step = 0
img_num = 35
 
dr = 4
radius = 370
 
num_theta= 2048
num_rad = 2*dr

num_var_t = 6
num_var_r = 5

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

#A0ft_svd_list = ringModel.generate_basis_matrices()

ringModel.l1_ratio = 1
ringModel.max_iters = 500
ringModel.compute_lipschitz(A0_stack)

synth_image_list = []
for j in range(3):
    synth_image = np.zeros((num_rad,num_theta))
    for i in range(5):
        rand_t    =   random.randint(0,num_var_t-1)
        rand_r    =   random.randint(0,num_var_r-1)
        rand_theta  = random.randint(0,num_theta-1)
        rand_rad    =   random.randint(0,num_rad-1)

        synth_image += 10*random.random()*EM.gaussian_basis_wrap_2D(num_theta,dtheta, rand_theta,  var_theta[rand_t], 
                                                           num_rad,  drad,   rand_rad,    var_rad[rand_r]).reshape((num_rad,num_theta))
    synth_image_list.append(synth_image)


partial_fit = partial(RM.fit_circulant_FISTA_Multiprocess,ringModel=ringModel,A0_stack=A0_stack,positive=1,benchmark=1,verbose=1)

rm_test = partial_fit(synth_image_list[0])
print(rm_test.times)

pool = multiprocessing.Pool(processes=3)
ringModels_synth = pool.map(partial_fit,synth_image_list)

for i in range(len(ringModels_synth)):
    print(ringModels_synth[i].times)

#np.save('ringModels_synth.npy',ringModels_synth)

