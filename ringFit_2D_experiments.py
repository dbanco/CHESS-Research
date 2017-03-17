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
import RingModel as RM
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

# Data, Interpolation, Fitting Parameters
load_step = 0
img_num = 35

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
    		              num_theta, num_rad, var_theta, var_rad, dtheta, drad)

# Compute basis matrix and interpolate ring to polar coordinates 
A0_stack = ringModel.generate_basis_matrix_stack()
polar_image = ringModel.compute_polar_image(sample)

ringModel.l1_ratio = 1
ringModel.max_iters = 500
ringModel.compute_lipschitz(A0_stack)

# Fit data
ringModel.fit_circulant_FISTA(A0_stack,positive=1,benchmark=1,verbose=1)



