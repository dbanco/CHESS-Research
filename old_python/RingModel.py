"""
Created on Thu Mar 16 04:12:40 2017

@author: Daniel Banco
"""

import numpy as np
from numpy.linalg import norm 

import DataAnalysis as DA
import EllipticModels as EM
import LassoSolvers
import RingImageProcessing as RingIP
import CirculantOperations as CO
import os
import time

class RingModel:
    def __init__(self, load_step, img_num, ring_radius_low, ring_radius_high, 
				   num_theta, num_rad, dtheta, drad, var_theta, var_rad):

        # Metadata                 
        self.load_step = load_step
        self.img_num = img_num        
        self.ring_radius_low = ring_radius_low
        self.ring_radius_high = ring_radius_high
        
        # Domain
        self.num_theta = num_theta
        self.num_rad = num_rad
        self.dtheta = dtheta
        self.drad = drad         
        self.var_theta = var_theta
        self.var_rad = var_rad


    def compute_polar_image(self, sample):
        image = sample.load_image(self.img_num, self.load_step)
		
        polar_image = np.zeros((self.num_rad, self.num_theta))

        for i, r in enumerate(np.arange(self.ring_radius_low,self.ring_radius_high,1)):
            fi, eta_domain = RingIP.azimuthal_projection(image,
												                sample.true_center,
                                                         r, 
                                                         0, 
                                                         2*np.pi-self.dtheta,
                                                         self.num_theta)
            polar_image[i,:] = fi

        self.polar_image = polar_image

    def generate_A0ft_svd_list(self):
        A0ft_list = EM.unshifted_basis_ft_svd_list(self.var_theta,
									                               self.var_rad,
									                               self.dtheta,
									                               self.drad,
									                               self.num_theta,
									                               self.num_rad)
        return A0ft_list

    def compute_lipschitz(self,A0_stack):
        lipschitz = 0
        for i in range(self.var_theta.shape[0]):
            for j in range(self.var_rad.shape[0]):
                lipschitz += max(np.linalg.eig( A0_stack[:,:,i,j] ).real)
        print(lipschitz)
        print(lipschitz*self.num_rad*self.num_theta)
        self.lipschitz = lipschitz

    def fit_circulant_FISTA(self,A0ft_stack,out_path,positive=1,benchmark=0,verbose=0):
        x_hat, times = LassoSolvers.fista_circulant_2D(A0ft_stack, self.polar_image, self.lipschitz, self.l1_ratio, self.max_iters, 
                         positive=positive,
											   benchmark=benchmark,
											   verbose=verbose) 

        y_hat = CO.Ax_ft_2D(A0ft_stack,x_hat)

        self.fit_image = y_hat  
        self.times = times      
        self.coefs = x_hat
        self.fit_error = norm(y_hat-self.polar_image)
        self.rel_fit_error = self.fit_error/norm(self.polar_image)
        
        np.save(os.path.join(out_path,'_load_',str(self.load_step),'_img_',str(self.img_num),'.npy'),self)

        
        return self

    def fit_circulant_svd_FISTA(self,A0ft_svd_list,positive=1,benchmark=0,verbose=0):
        x_hat, times = LassoSolvers.fista_circulant_2D_SVD(A0ft_svd_list,len(self.var_theta), len(self.var_rad), self.polar_image,self.lipschitz, self.l1_ratio, self.max_iters, 
                         positive=positive,
											   benchmark=benchmark,
											   verbose=verbose) 

        y_hat = LassoSolvers.Ax_ft_2D_SVD(A0ft_svd_list,x_hat)

        self.fit_image = y_hat  
        self.times = times      
        self.coefs = x_hat
        self.fit_error = norm(y_hat-self.polar_image)
        self.rel_fit_error = self.fit_error/norm(self.polar_image)

    def fit_parallel_circulant_FISTA(self,A0ft_svd_list,positive=1,benchmark=0,verbose=0):
        x_hat, times = LassoSolvers.fista_circulant_2D_Parallel_SVD(A0ft_svd_list,len(self.var_theta), len(self.var_rad), self.polar_image,
 self.lipschitz, self.l1_ratio, self.max_iters,   positive=positive,

           benchmark=benchmark,

           verbose=verbose)

        y_hat = CO.Ax_ft_2D_Parallel_SVD(A0ft_svd_list,x_hat)

        self.fit_image = y_hat
        self.times = times
        self.coefs = x_hat
        self.fit_error = norm(y_hat-self.polar_image)
        self.rel_fit_error = self.fit_error/norm(self.polar_image)


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
    
def compute_lipschitz(ringModel,A0_stack):
    lipschitz = 0
    for i in range(ringModel.var_theta.shape[0]):
        for j in range(ringModel.var_rad.shape[0]):
            eigs,vecs = np.linalg.eig( np.dot(A0_stack[:,:,i,j],A0_stack[:,:,i,j].T) )
            lipschitz += np.max(eigs.real)
    print(lipschitz)
    print(lipschitz*ringModel.num_rad*ringModel.num_theta)
    ringModel.lipschitz = lipschitz
    
def generate_A0_stack(ringModel):
    A0_stack = EM.unshifted_basis_matrix_stack(ringModel.var_theta,
							                               ringModel.var_rad,
							                               ringModel.dtheta,
							                               ringModel.drad,
							                               ringModel.num_theta,
							                               ringModel.num_rad)
    return A0_stack

def fit_circulant_FISTA_Multiprocess(ringModel,A0ft_stack,out_path,positive=1,benchmark=0,verbose=0):
    x_hat, times = LassoSolvers.fista_circulant_2D(A0ft_stack, ringModel.polar_image, ringModel.lipschitz, ringModel.l1_ratio, ringModel.max_iters, 
                     positive=positive,
									   benchmark=benchmark,
									   verbose=verbose) 

    y_hat = CO.Ax_ft_2D(A0ft_stack,x_hat)

    ringModel.fit_image = y_hat  
    ringModel.times = times      
    ringModel.coefs = x_hat
    ringModel.fit_error = norm(y_hat-ringModel.polar_image)
    ringModel.rel_fit_error = ringModel.fit_error/norm(ringModel.polar_image)

    np.save(out_path+'_load_'+str(ringModel.load_step)+'_img_'+str(ringModel.img_num)+'.npy',ringModel)
