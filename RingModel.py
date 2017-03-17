"""
Created on Thu Mar 16 04:12:40 2017

@author: Daniel Banco
"""

import numpy as np
import matplotlib.pyplot as plt

from numpy.linalg import norm 

import DataAnalysis as DA
import EllipticModels as EM
import LassoSolvers as Lasso
import os
import time

class RingModel:
    def __init__(self, load_step, img_num, ring_radius_low, ring_radius_high, 
				 polar_image, num_theta, num_rad, dtheta, drad, var_theta, var_rad):

        # Metadata                 
        self.load_step = load_step
        self.img_num = img_num        
        self.ring_radius_low = ring_radius
		self.ring_radius_high = ring_radius_high
        
        # Interpolated data
        self.polar_image = polar_image
        
		# Domain
		self.num_theta = num_theta
		self.num_rad = num_rad
		self.dtheta = dtheta
		self.drad = drad         
		self.var_theta = var_theta
        self.var_rad = var_rad
      
		self.l1_ratio = 1
		self.lipschitz = 1
		self.max_iters = 500

	def compute_polar_image(self, sample):
    	image = specimen.load_image(self.img_num, self.load_step)
		
		polar_image = np.zeros((self.num_rad, self.num_theta))

		for i, r in enumerate(np.arange(self.ring_radius_low,self.ring_radius_high,1)):
			fi, theta_domain = RingIP.azimuthal_projection(image,
												           sample.true_center,
														   r, 0, 2*np.pi-self.dtheta,
                                                           self.num_theta)
			polar_image[i,:] = fi

		self.polar_image = polar_image

	def generate_basis_matrix_stack(self):
		A0_stack = generate_basis_matrix_stack(self.var_theta,
									           self.var_rad,
									           self.dtheta,
									           self.drad,
									           self.num_theta,
									           self.num_rad)
		return A0_stack

	def compute_lipschitz(self,A0_stack):
        flat_shape = (self.num_rad*self.var_theta.shape[0]*self.var_rad.shape[0],self.num_theta)
		eig = np.linalg.eig( np.dot(A0_stack.reshape(flat_shape).T,
                                    A0_stack.reshape(flat_shape) ))
		lipschitz = eig[0].real*num_rad*num_theta/500
		self.lipschitz = lipschitz


	def fit_circulant_FISTA(self,A0_stack,positive=1,benchmark=0,verbose=0)
		x_hat, times = LassoSolvers.fista_circulant_2D_Parallel(B0_stack, self.polar_image, 
				                               self.lipschitz, self.l1_ratio, self.max_iters, 
				                               positive=positive,
											   benchmark=benchmark,
											   verbose=verbose) 

		y_hat = Lasso.Ax_ft_2D(A0_stack,x_hat)

		self.fit_image = y_hat  
		self.times = times      
		self.coefs = x_hat
        self.fit_error = norm(y_hat-polar_image)
        self.rel_fit_error = self.fit_error/norm(polar_image)

	def print_fit_stats(self):
		
		print('Relative Error: '  + str(rel_fit_error))

		pos_coef = np.sum(self.coefs>0)
		tot_coef = len(self.coefs.ravel())
		sparsity0 = pos_coef/tot_coef
		one_coef = np.sum(self.coefs>1)
		sparsity1 = one_coef/tot_coef
		
		print('x > 0         |  x > 1             | total x')
		print(str(pos_coef) + '       |  ' + str(one_coef) + '             | ' + str(tot_coef))
		print(str(sparsity0) + ' |  ' + str(sparsity1))
    
