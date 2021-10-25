## -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 13:35:59 2017

@author: dbanco02
"""

## Init
from __future__ import division
from __future__ import print_function

import numpy.distutils.system_info as sysinfo
sysinfo.get_info('atlas')

import numpy as np
import os
import RingModel as RM
import EllipticModels as EM
import LassoSolvers
import CirculantOperations as CO
import random
from functools import partial
import multiprocessing
import time

# Data, Interpolation, Fitting Parameters
load_step = 0
img_num = 35
 
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
A0_stack = EM.unshifted_basis_matrix_stack(var_theta,var_rad,
                                           dtheta,
                                           drad,
                                           num_theta, 
                                           num_rad)
#%% Test individual operations
reload(CO)
a = time.time()
np.fft.fft2(A0_stack[:,:,0,0])

#CO.AtR_ft_2D(A0_stack, A0_stack[:,:,0,0])
print(time.time()-a)


#%% Test Algorithm
a = time.time()
CO.Ax_ft_2D(A0_stack, A0_stack)
#LassoSolvers.fista_circulant_2D(A0_stack, np.sum(np.sum(A0#_stack,2),2), 6e5, 1, 5, eps=10**(-7), positive=1, verbose=1, #benchmark=1)
print(time.time()-a)

#np.save('ringModels_synth.npy',ringModels_synth)

