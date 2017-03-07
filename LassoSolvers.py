# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 16:13:12 2017

Library of functions for solving Lasso problem

@author: Daniel Banco
"""

## Init
from __future__ import division
from __future__ import print_function

import numpy as np
import os
import time

from sklearn.linear_model import Lasso
from numpy.linalg import norm 
import matplotlib.pyplot as plt
from scipy.signal import argrelmax
from math import sqrt

import RingImageProcessing as RingIP
import EllipticModels as EM

"Soft thresholding function with vector input and output"
def soft_thresh(x, l):
    return np.sign(x) * np.maximum(np.abs(x) - l, 0.)    

"FISTA algorithm"
def fista(A, b, l, maxit, positive=0, benchmark=0):
    
    if benchmark: 
        start_time = time.time()
        
    x = np.zeros(A.shape[1])
    t = 1
    z = x.copy()
    L = norm(A)**2
    for it in range(maxit):
        xold = x.copy()
        
        # Project onto constraint
        z = z + A.T.dot(b - A.dot(z))/L
        # Enforce positivity 
        if positive:
            z[z<0] = 0
        # Sot Threshold
        x = soft_thresh(z,1/L)
        # Update Step size
        t0 = t
        t = (1. + sqrt(1. + 4. * t ** 2)) / 2.
        z = x + ((t0 - 1.) / t) * (x - xold)
        if positive:
            z[z<0] = 0
             
    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x

def LASSO_cost(A,b,x):
    """
    Compute LASSO cost function
    """
    return 0.5/len(b)*np.sum((A.dot(x) - b)**2)
    
def LASSO_grad(A,b,x):
    """
    Compute gradient of LASSO cost function
    """
    return A.T.dot(A.dot(x)-b)/len(b)


def prox_ADMM(a, lam, max_iters, eps):
    """
    Solves
    argmin_x 1/2||x-a||^2_2 + lambda*(||x||_1 + Ind(x)_+)    
    
    Inputs:
    v           Vector
    lam          Soft thresholding parameter
    max_iters   Max number of iterations
    delta       Relative signal change after NPG iteration
    eta         ADMM convergence tuning constant
    
    Outputs:
    alpha       Solution to minimization problem (x)
    """
        
    # Init
    alpha = a.copy()
    alpha[a<0] = 0
    z = np.zeros(a.shape)    
    v = a - alpha

    # First iteration
    z_new = soft_thresh(alpha - v,lam)
    alpha_new = 0.5*(a + z_new + v)
    alpha_new[alpha_new<0] = 0
    v = v + z_new - alpha_new        
    
    # Update 
    alpha_old = alpha.copy()
    alpha = alpha_new.copy()
    z_old = z.copy()
    z = z_new.copy()

    #rho = 1
    for j in range(max_iters):
        
        """
        a_bar = (a + rho*(s_1 + v_1))/(1+rho)
        alpha = 1
        """
        
        # Convergence check
        conv_1 = norm(z - z_old)/norm(z)
        conv_2 = norm(alpha- alpha_old) /norm(alpha)  
        if max(conv_1,conv_2) < eps:
            break
        
        # Iterate
        z_new = soft_thresh(alpha - v,lam)
        alpha_new = 0.5*(a + z_new + v)
        alpha_new[alpha_new<0] = 0
        v = v + z_new - alpha_new          
        
        # Update
        alpha_old = alpha.copy()
        alpha = alpha_new.copy()
        z_old = z.copy()
        z = z_new.copy()
        
    print(j)
    return alpha
 
def prox_LASSO(A,b,a, lam, L, max_iters, eps):
    """ 
    Solves
    argmin_x 1/2||x-a||^2_2 + lambda*(||x||_1 + Ind(x)_+)    
    
    Inputs:
    v           Vector
    lam         Soft thresholding parameter
    max_iters   Max number of iterations
    delta       Relative signal change after NPG iteration
    eta         ADMM convergence tuning constant
    
    Outputs:
    alpha       Solution to minimization problem (x)
    """
        
    # Init
    alpha = a.copy()
    alpha[a<0] = 0

    for j in range(max_iters):
                
        # Iterate
        alpha_new = alpha - LASSO_grad(A,b,alpha)/L
        alpha_new[alpha_new<0] = 0        
        
        # Update
        alpha = alpha_new.copy()

    return alpha 
 
def NN_LASSO_obj(A,x,b,lam):
    """
    Evaluate nonnegative LASSO objective function
    ||Ax-b||^2_2 + lam||x||_1 + Ind(x)
    """
    
    if np.sum(x<0):
        print('Negative coefficient found')
        return float("inf")
    else:
        return norm(A.dot(x)-b)**2 + lam*np.sum(np.abs(x))

def NN_LASSO_obj_quad_approx(L,lam,x,y,A,b):
    """
    Evaluate nonnegative LASSO quadratic approximation
    f(y) + <x-y,gradf(y)> + 0.5*L*norm(x-y)**2 + g(x)
    """
    if np.sum(x<0):
        print('Negative coefficient found')
        return float("inf")
    else:
        Q_L = LASSO_cost(A,b,y) + \
              np.dot(x-y,LASSO_grad(A,b,y)) + \
              0.5*L*norm(x-y)**2 + lam*np.sum(np.abs(x))
        return Q_L


"""
def Quadratic_approx(beta,x,x_bar,A,b):
    ""
    Evaluate nonnegative LASSO quadratic approximation
    f(x_bar) + <x-x_bar,gradf(y)> + 0.5/beta*norm(x-x_bar)**2
    ""
    if np.sum(x<0):
        print('Negative coefficient found')
        return float("inf")
    else:
        Q_L = LASSO_cost(A,b,x_bar) + \
              np.dot(x-x_bar,LASSO_grad(A,b,x_bar)) + \
              0.5/beta*norm(x-x_bar)**2
        return Q_L

def PNPG(x0,u,gamma,b,n,m,xc,eta,eps,max_iters):
    
    x_1 = np.zeros(x0.shape)
    i = 0
    kappa = 0
    
    # Initialize by BB method
    grad1 = LASSO_grad(A,y,x0)
    temp = 1e-5*grad1/norm(grad1)
    grad2 = LASSO_grad(A,y,x0-temp)
    beta = np.abs(np.dot(grad1-grad2,temp))/np.sum((opt-x0)**2)
    
    for it in range(max_iters):
        i += 1
        kappa += 1
        x_2 = x_1.copy()
        
        # Backtracking search
        while True:
            B = beta_old/beta
            if(i > 1):
                theta = 1/gamma + np.sqrt(b + B*theta_1**2)
            else:
                theta = 1
            T_const = (theta_1 - 1)/theta
            x_bar = max(x_1 + T_const*(x_1 - x_2),0)
            x = prox(x_bar - beta*LASSO_grad(A,y,x_bar))
            

"""

#def linearModel(x,A,y,b):
#    Ax = A.dot(x)
#
#    g = Ax - b
#    f = np.sum(g**2)/2
#    h = np.ones(x.shape)
#    
#    return f, A.T.dot(g), hessian1(x,A,opt,h)
    
def hessian1(x,A,b,opt=2):
    h = np.ones(x.shape)
    yy = A.dot(x)
    if(opt == 1):
        hh = A.T.dot(h*yy)
    else:
        hh = np.dot(yy,h*yy) 
    return hh

#def nonegPenalty(x):
#    temp = x < 0
#    f = np.dot(x(temp),x(temp))
#    
#    g = np.zeros(x.shape)
#    g(temp) = 2*x(temp)
#    
#    hh = hessian2(x,A,opt)
    
def hessian2(x,opt):
    temp = x < 0
    if(opt == 1):
        hh = np.zeros(x.shape)
        hh[temp,:] = x[temp,:]*2
    else:
        x = np.reshape(x,len(temp[:]),np.array([]))
        y = x[temp,:]
        hh = np.dot(y,y*2)
    return(hh)

def hessian(x,opt):
    hh = 0
    for j in range(2):
        hh = hh + hessian1(x,A,bopt)*coef1 + hessian2(x,opt)*coef2

"NN FISTA algorithm"
def nn_fista(x0, A, b, lam, L, maxit, eps, benchmark=0):
    
    if benchmark: 
        start_time = time.time()
        
    # Step 0 
    stepShrink = 0.5
    preAlpha = 0
    thresh = 1e-4
    theta = 0
    cumu = 0
    cumuTol = 4
    oldCost = 0
    x = x0
    for it in range(maxit):
        temp = (1. + sqrt(1. + 4. * t ** 2)) / 2. 
        y = x + ((theta - 1.) / theta) * (x - x_old)   
        
        theta = temp
        x_old = x
        
        # Initialize theta step size if not given #############################
        if(theta==-1):
            oldCost, grad, hessian = obj_func_hess(x,A,b)
            t = hessian(grad,2)/np.dot(grad,grad)
        else:
            oldCost, grad = obj_func(x,A,b)
        ############################################################
        
        # Begin line search
        i = 0
        while True:
            newX = y - grad/t
            newX[newX>0] = 0
            newCost = obj_func(newX)[0]
            if(newCost<=oldCost + grad.dot(newX-y) + theta/2*norm(newX-y)**2):
                break
            else:
                t = t/stepShrink
        if(i == 1):
            cumu = cumu + 1
            if(cumu >= cumuTol):
                t = t*stepShrink
                cumu = 0
        else:
            cumu = 0
            
        difx = norm(newX-x)/norm(newX)
        x = newX
        cost = newCost
        if(difx<thresh):
            break        
    
    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x


import scipy.optimize as opt
def obj_func(x,A,b,lam,eps):
    """
    Computes full LASSO objective function and Jacobian
    """
    
    obj = LASSO_cost(A,b,x)
    obj += lam*np.sum(np.sqrt(x**2 + eps))
    
    grad = A.T.dot(A.dot(x)-b)/len(b)
    grad += lam*2*x/np.sqrt(x**2 + eps)
    
    return obj, grad

def obj_func_circ(x,A0_ft,A0_tran_ft,b,maxima,num_var,m,n,lam,eps):
    """
    Compute LASSO cost function leveraging circulant structure
    """
    res = Ax_ft(A0_ft,x,maxima,m,n) - b
    obj = 0.5/len(b)*np.sum(res**2) + lam*np.sum(np.abs(x))   
    grad = ATb_ft(A0_tran_ft,res,num_var,maxima)/len(b) + lam*2*x/np.sqrt(x**2 + eps)
    
    return obj, grad

def obj_func_circ_callback(x,A0_ft,A0_tran_ft,b,maxima,num_var,m,n,lam,eps):
    """
    Compute LASSO cost function leveraging circulant structure
    """
    res = Ax_ft(A0_ft,x,maxima,m,n) - b
    obj = 0.5/len(res)*np.sum(res**2) + lam*np.sum(np.abs(x)) 
    obj_hist = np.load('obj_hist.npy')
    obj_hist = np.append(obj_hist,obj)
    np.save('obj_hist.npy',obj_hist)

def LASSO_approx_tnc(A0,A0_tran,maxima,num_var,num_theta,b,lam,eps,max_iter,callback=0,benchmark=0):
  
    if benchmark: 
        start_time = time.time()
    
    x0 = np.zeros(num_var*len(maxima))
    m = num_theta
    n = A0.shape[1]
    A0_ft = np.fft.fft(A0,axis=0)    
    A0_tran_ft = np.fft.fft(A0_tran,axis=0)  
    
    # Construct positivity bounds
    bnd = []
    for i in range(len(x0)):
        bnd.append((0,None))
        
    if callback:
        obj_hist = np.array(())
        np.save('obj_hist.npy',obj_hist)
        callbackf = lambda x: obj_func_circ_callback(x,A0_ft,A0_tran_ft,
                                                     b,maxima,num_var,m,n,
                                                     lam,eps)
        
        x = opt.fmin_l_bfgs_b(obj_func_circ,
                                   x0,
                                   args=(A0_ft,A0_tran_ft,b,maxima,num_var,m,n,lam,eps),
                                   bounds=bnd,
                                   maxiter=max_iter,
                                   callback=callbackf)
    else:
        x = opt.fmin_l_bfgs_b(obj_func_circ,
                       x0,
                       args=(A0_ft,A0_tran_ft,b,maxima,num_var,m,n,lam,eps),
                       bounds=bnd,
                       maxiter=max_iter) 
    
    if benchmark:
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x


"Vector formatting helper function for circulant matrix-vector product"
def x_to_x_ft(x, maxima, m, n):
    """ 
    Inputs:
        x       coefficient vector
        maxima  vector of locations of maxima in signal
        m       number of rows in x_ft
        n       number of columns in x_ft
    Output:
        x_ft    Each column is the fourier transform of the coefficients
                corresponding to the circulant matrix Ai
    """
    x_pad = np.zeros((m,n))
    x_ft = np.zeros((m,n))
    num_maxima = len(maxima)
    for i in range(n):
        x_pad[maxima,i] = x[i*num_maxima:(i+1)*num_maxima]      
        
    x_ft = np.fft.fft(x_pad,axis=0)
    
    return x_ft

"Circulant matrix-vector product subroutine"
def Ax_ft(A0ft, x, maxima, m, n):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        x       Coefficient vector
        maxima  vector of locations of maxima in signal
        m       number of rows in x_ft
        n       number of columns in x_ft
    """
    Ax = np.zeros((A0ft.shape[0]))
    
    x_ft = x_to_x_ft(x,maxima,m,n)
    
    for ii in range(A0ft.shape[1]):
        Ax += np.fft.ifft(np.multiply(A0ft[:,ii],x_ft[:,ii])).real
        
    return Ax

"Transposed circulant matrix-vector product subroutine"
def ATb_ft(A0_tran_ft,R,num_var,maxima):
    # Inputs:
    #    A0_tran_ft    Each column is the first row of a circulant matrix Ai
    #    R             Residual vector 
    # Output:
    #    ATb           Remember that this will be zero padded

    num_maxima = len(maxima)
    ATb = np.zeros((num_var*num_maxima))
    R_ft = np.fft.fft(R)          
    
    # num_vars    
    for ii in range(num_var):
        ATb[ii*num_maxima:(ii+1)*num_maxima] = np.fft.ifft(np.multiply(A0_tran_ft[:,ii],R_ft)).real[maxima]   

    return ATb
    
"FISTA algorithm using circulant matrix-vector product subroutines"
def fista_circulant(A, A0, A0_tran, maxima, b,
                    l1_ratio, maxit, positive=0, benchmark=0):
    
    if benchmark: 
        start_time = time.time()
            
    m = A.shape[0]    
    n = A0.shape[1]
    A0ft = np.fft.fft(A0,axis=0)    
    A0_tran_ft = np.fft.fft(A0_tran,axis=0)  
    
    x = np.zeros(A.shape[1])
    t = 1
    z = x.copy()
    L = norm(A)**2
    for it in range(maxit):
        xold = x.copy()
        
        # Arrange x coefficents as matrix in fourier domain 
        R = b - Ax_ft(A0ft, z, maxima, m, n)
        z = z + ATb_ft(A0_tran_ft,R,n,maxima)/L
        
        # Enforce positivity on coefficients
        if positive:
            z[z<0]=0

        x = soft_thresh(z,1/L)
        
        t0 = t
        t = (1. + sqrt(1. + 4. * t ** 2)) / 2.
        z = x + ((t0 - 1.) / t) * (x - xold)
        
        if positive:
            z[z<0]=0
    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x
    
"Circulant matrix-vector product subroutine"
def Ax_ft_2D(A0ft, x):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        x       Coefficient vector
    """
    Ax = np.zeros((A0ft.shape[0],A0ft.shape[1]))
    
    x_ft = np.fft.fft2(x,axes=(0,1))
    
    for tv in range(A0ft.shape[2]):
        for rv in range(A0ft.shape[3]):
            Ax += np.fft.ifft2(A0ft[:,:,tv,rv]*x_ft[:,:,tv,rv]).real
        
    return Ax

from multiprocessing import Pool

"Circulant matrix-vector product subroutine"
def AtR_ft_2D_Parallel(A0ft, R, P):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        x       Coefficient vector
        P       Number of processes to run simultaneously
    """
    AtR = np.zeros((A0ft.shape))
    
    R_ft = np.fft.fft2(R)

    pool = Pool(processes=P)
    
    n = A0ft.shape[2]
    m = A0ft.shape[3]
    k = 0
    index = np.zeros((n*m,2))
    for tv in range(n):
        for rv in range(m):
            index[k,0] = tv
            index[k,1] = rv
            k += 1
    pool.map        
    AtR[:,:,tv,rv] += np.fft.ifft2(A0ft[:,:,tv,rv]*R_ft).real
        
    return AtR
    

"Circulant matrix-vector product subroutine"
def Ax_ft_2D_Parallel(A0ft, x, P):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        x       Coefficient vector
        P       Number of processes to run simultaneously
    """
    Ax = np.zeros((A0ft.shape[0],A0ft.shape[1]))
    
    x_ft = np.fft.fft2(x,axes=(0,1))
    
    for tv in range(A0ft.shape[2]):
        for rv in range(A0ft.shape[3]):
            Ax += np.fft.ifft2(A0ft[:,:,tv,rv]*x_ft[:,:,tv,rv]).real
        
    return Ax

"Circulant matrix-vector product subroutine"
def AtR_ft_2D(A0ft, R):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        x       Coefficient vector
    """
    AtR = np.zeros((A0ft.shape))
    
    R_ft = np.fft.fft2(R)

    for tv in range(A0ft.shape[2]):
        for rv in range(A0ft.shape[3]):
            AtR[:,:,tv,rv] += np.fft.ifft2(A0ft[:,:,tv,rv]*R_ft).real
        
    return AtR


"FISTA algorithm using circulant matrix-vector product subroutines"
def fista_circulant_2D(A0, b, L, l1_ratio, maxit, eps=10**(-8), positive=0, verbose=0, benchmark=0,):
    # A0 is a bunch of slices indexed by variance and radius
    if benchmark: 
        start_time = time.time()
            
        
#    num_rad = A0.shape[0]
#    num_theta = A0.shape[1]
    num_var_theta = A0.shape[2]
    num_var_rad = A0.shape[3]
    Linv = 1/L
    # Compute fft of each slice
    A0ft = np.zeros(A0.shape)
    for tv in range(num_var_theta):
        for rv in range(num_var_rad):
            A0ft[:,:,tv,rv] = np.fft.fft2(A0[:,:,tv,rv]) 
            
    
    x = np.zeros(A0.shape)
    t = 1
    z = x.copy()

    for it in range(maxit):
        xold = x.copy()
        
        # Arrange x coefficents as matrix in fourier domain 
        R = b - Ax_ft_2D(A0ft,z)
        z = z + AtR_ft_2D(A0ft,R)*Linv
        
        # Enforce positivity on coefficients
        if positive:
            z[z<0]=0

        x = soft_thresh(z,l1_ratio*Linv)
        
        t0 = t
        t = (1. + sqrt(1. + 4. * t ** 2)) / 2.
        z = x + ((t0 - 1.) / t) * (x - xold)
        
        if positive:
            z[z<0]=0
        
        criterion = np.sum(np.abs(x - xold))/len(x.ravel())
        
        if verbose:
            res = np.sum( (b.ravel() - Ax_ft_2D(A0ft,x).ravel())**2 )
            L1 =  l1_ratio*np.sum(np.abs( x.ravel() ))
            obj =  res + L1
            print('Iteration  ' +\
                  'Objective      ' +\
                  'Relative Error      ' +\
                  'Residual            ' +\
                  'L1 Penalty          ' +\
                  'Criterion' )
            
            print( str(it) +' of '+ str(maxit) + '   ' +\
                   str(obj)                    + ' ' +\
                   str(np.sqrt(res)/norm(b))   + '       ' +\
                   str(res)                    + '       ' +\
                   
                   str(L1)                     + '       ' +\
                   str(criterion))
#            print('Iter: '     + str(it) +' of '+ str(maxit) +\
#                  ', Obj: '    + str(obj) +\
#                  ', Res: '    + str(res) +\
#                  ', RelErr: ' + str(np.sqrt(res)/norm(b))
#                  ', L1: '     + str(L1)   +\
#                  ', Crit: '   +str(criterion))
        if(criterion < eps):
            break
            
        
    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x
    
def coord_ascent(A,
                 A_tran,
                 A0,
                 A0_tran,
                 maxima,
                 b,
                 alpha, 
                 max_iter,
                 positive=0, 
                 benchmark=0):
    
    if benchmark: 
        start_time = time.time()
            
    m = A.shape[0]
    n = A.shape[1]
    num_vars = A0.shape[1]
    num_maxima = maxima.shape[0]

    # precompute sum of sqaures of columns
    z = np.sum(A**2,0) 
   
    # initialize x
    x = np.zeros(n)
    x_old = np.zeros(num_maxima)             
    AtR = np.zeros(m)
    R = np.zeros(m) 


    rho = np.zeros(num_maxima)   
    tol = 1e-4
    gap = tol + 1.0
    d_x_tol = tol       
       
    # R = b - np.dot(A,x) but x is zeros   
    R = b.copy()
    
    # tol *= np.dot(b,b)
    tol *= np.dot(b,b)
 
    for iters in range(max_iter):
        x_max = 0.0
        d_x_max = 0.0
        
        for ii in range(num_vars):
            # 1)            
            x_old = x[ii*num_maxima:(ii+1)*num_maxima]
            
            # Remove contribution to residual of current coefficient
            # R += x_ii*A[:,ii]
            R += np.dot(A[:,ii*num_maxima:(ii+1)*num_maxima],x[ii*num_maxima:(ii+1)*num_maxima])
#            plt.figure(1)
#            plt.plot(R)
#            plt.pause(0.25)
            # 2)
            # Soft threshold update with optional positivity constraint
            rho = np.dot(A_tran[ii*num_maxima:(ii+1)*num_maxima,:],R)        
#            plt.figure(2)
#            plt.plot(rho,'o')
#            plt.pause(0.25)            
            # rho = np.dot(A[:,ii],R)
            if positive:
                rho[rho<0]=0
            x[num_maxima*ii:num_maxima*(ii+1)] = np.sign(rho)*np.maximum(np.abs(rho)-alpha,np.zeros((len(rho))))/z[ii]
            
            # 3)            
            # Add contribution to residual of updated coefficient
            # R -= x[ii]*A[:,ii]
            R -= np.dot(A[:,ii*num_maxima:(ii+1)*num_maxima],x[ii*num_maxima:(ii+1)*num_maxima])
#            plt.figure(1)            
#            plt.plot(R)
#            plt.pause(0.25)
            # update the maximum absolute coefficient update
            d_x_ii = np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima] - x_old))
            if d_x_ii > d_x_max:
                d_x_max = d_x_ii

            if np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima])) > x_max:
                x_max = np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima]))

#            x_pad = np.zeros(m)
#            x_pad[maxima] = x[ii*num_maxima:(ii+1)*num_maxima]
#            b_hat = np.zeros(m)
#            for i in range(num_vars):
#                b_hat += np.fft.ifft(np.multiply(Aft[:,i],np.fft.fft(x_pad))).real
#            
#            plt.figure(3)
#            plt.plot(b_hat)
#            plt.pause(0.25)
#        
        if (x_max == 0.0 or
            d_x_max / x_max < d_x_tol or
            iters == max_iter - 1):
            # the biggest coordinate update of this iteration was smaller
            # than the tolerance: check the duality gap as ultimate
            # stopping criterion

            # AtR = np.dot(A.T, R)
            AtR = np.dot(A_tran, R)


            if positive:
                dual_norm_AtR = np.max(AtR)
            else:
                dual_norm_AtR = np.max(np.abs(AtR))

            R_norm2 = np.dot(R, R)

            x_norm2 = np.dot(x, x)

            if (dual_norm_AtR > alpha):
                const = alpha / dual_norm_AtR
                A_norm2 = R_norm2 * (const ** 2)
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            l1_norm = np.sum(x)

            # np.dot(R.T, b)
            gap += (alpha * l1_norm - const * np.dot(R.T,b))

            if gap < tol:
                break
                # return if we reached desired tolerance
                
    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x

def coord_ascent_circ(A,
                      A_tran,
                      A0,
                      A0_tran,
                      maxima,
                      b,
                      alpha, 
                      max_iter,
                      positive=0, 
                      benchmark=0):

    if benchmark: 
        start_time = time.time()

    m = A.shape[0]
    n = A.shape[1]
    num_vars = A0.shape[1]
    num_maxima = maxima.shape[0]

    # precompute sum of sqaures of columns
    z = np.sum(A**2,0) 
   
    # initialize x
    x = np.zeros(n)
    x_old = np.zeros(num_maxima)       
    x_pad = np.zeros(m)        
    AtR = np.zeros(m)
    R = np.zeros(m) 
    Rft = np.zeros(m) 


    rho = np.zeros(num_maxima)   
    tol = 1e-4
    gap = tol + 1.0
    d_x_tol = tol  
    

    Aft = np.fft.fft(A0,axis=0)
    A0_tran_ft = np.fft.fft(A0_tran,axis=0)       
       
    # R = b - np.dot(A,x) but x is zeros   
    R = b.copy()
    Rft = b.copy()
    
    # tol *= np.dot(b,b)
    tol *= np.dot(b,b)
 
    for iters in range(max_iter):
        x_max = 0.0
        d_x_max = 0.0
        
        for ii in range(num_vars):
            # 1)            
            x_old = x[ii*num_maxima:(ii+1)*num_maxima]
            x_pad = np.zeros(m)
            x_pad[maxima] = x_old
            x_ft = np.fft.fft(x_pad)
            
            # Remove contribution to residual of current coefficient
            # R += x_ii*A[:,ii]
            R += np.fft.ifft(np.multiply(Aft[:,ii],x_ft)).real
#            plt.figure(1)
#            plt.plot(R)
#            plt.pause(0.25)
            # 2)
            # Soft threshold update with optional positivity constraint
            rho = np.fft.ifft(np.multiply(A0_tran_ft[:,ii],np.fft.fft(R))).real[maxima]          
#            plt.figure(2)
#            plt.plot(rho,'o')
#            plt.pause(0.25)            
            # rho = np.dot(A[:,ii],R)
            if positive:
                rho[rho<0]=0
            x[num_maxima*ii:num_maxima*(ii+1)] = np.sign(rho)*np.maximum(np.abs(rho)-alpha,np.zeros((len(rho))))/z[ii]
            
            # 3)
            x_pad = np.zeros(m)
            x_pad[maxima] = x[ii*num_maxima:(ii+1)*num_maxima]
            x_ft = np.fft.fft(x_pad)
            
            # Add contribution to residual of updated coefficient
            # R -= x[ii]*A[:,ii]
            R -= np.fft.ifft(np.multiply(Aft[:,ii],x_ft)).real
#            plt.figure(1)            
#            plt.plot(R)
#            plt.pause(0.25)
            # update the maximum absolute coefficient update
            d_x_ii = np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima] - x_old))
            if d_x_ii > d_x_max:
                d_x_max = d_x_ii

            if np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima])) > x_max:
                x_max = np.max(np.abs(x[ii*num_maxima:(ii+1)*num_maxima]))

#            x_pad = np.zeros(m)
#            x_pad[maxima] = x[ii*num_maxima:(ii+1)*num_maxima]
#            b_hat = np.zeros(m)
#            for i in range(num_vars):
#                b_hat += np.fft.ifft(np.multiply(Aft[:,i],np.fft.fft(x_pad))).real
#            
#            plt.figure(3)
#            plt.plot(b_hat)
#            plt.pause(0.25)
#        
        if (x_max == 0.0 or
            d_x_max / x_max < d_x_tol or
            iters == max_iter - 1):
            # the biggest coordinate update of this iteration was smaller
            # than the tolerance: check the duality gap as ultimate
            # stopping criterion

            # AtR = np.dot(A.T, R)
            AtR = np.zeros(m)
            Rft = np.fft.fft(R)
            for i in range(num_vars):
                AtR += np.real(np.fft.ifft(np.multiply(A0_tran_ft[:,i],Rft)))


            if positive:
                dual_norm_AtR = np.max(AtR)
            else:
                dual_norm_AtR = np.max(np.abs(AtR))

            R_norm2 = np.dot(R, R)

            x_norm2 = np.dot(x, x)

            if (dual_norm_AtR > alpha):
                const = alpha / dual_norm_AtR
                A_norm2 = R_norm2 * (const ** 2)
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            l1_norm = np.sum(x)

            # np.dot(R.T, b)
            gap += (alpha * l1_norm - const * np.dot(R.T,b))

            if gap < tol:
                break
                # return if we reached desired tolerance
                
    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x
        
import circulant_lasso_ca as circ
def cython_circulant_coord_ascent(A,
                                  A0,
                                  maxima,
                                  b,
                                  alpha, 
                                  max_iter,
                                  positive=0, 
                                  benchmark=0):
    if benchmark: 
        start_time = time.time()
        
    num_theta = A0.shape[0]
    num_var   = A0.shape[1]    

    A_tran = np.transpose(A)

    # Create transposed convolution matrix
    A0_tran = np.zeros([num_theta,num_var])
    for i in range(num_var): 
        A0_tran[0,i] = A0[0,i]
        A0_tran[1::,i] = np.flipud(A0[1::,i])
    
    A0ft = np.fft.fft(A0,axis=0)
    A0_tranft = np.fft.fft(A0_tran,axis=0)    
    
    # Enforce Fortran contiguous storage
    A = np.array(A, order='F', copy=True)
    A_tran = np.array(A_tran, order='F', copy=True)
    A0 = np.array(A0, order='F',copy=True)
    A0_tran = np.array(A0_tran, order='F',copy=True)                                     
    A0ft = np.array(A0ft, order='F',copy=True)
    A0_tranft = np.array(A0_tranft, order='F',copy=True)

    # Call Cython function
    x = circ.lasso_solver(A, A0ft, A0_tran_ft, maxima, b, 
                           alpha, max_iter, positive=positive)
        
    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x

# This may be good for a sanity check, but not good for much else
import lasso_coord_ascent as coord      
def cython_coord_ascent(A,
                        b,
                        alpha, 
                        max_iter,
                        positive=0, 
                        benchmark=0):
    if benchmark:  
        start_time = time.time()
        
    # Enforce Fortran contiguous storage
    A = np.array(A, order='F', copy=True)

    # Call Cython function
    x = coord.lasso_solver(A, b, 
                           alpha, max_iter, positive)
        
    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x

#import circulant_fista as cf      
#def cython_circulant_fista(A,
#                           A0,
#                           maxima,
#                           b,
#                           alpha, 
#                           max_iter,
#                           positive=0, 
#                           benchmark=0):
#    if benchmark:  
#        start_time = time.time()
#        
#    num_theta = A0.shape[0]
#    num_var   = A0.shape[1]    
#
#    A_tran = np.transpose(A)
#
#    # Create transposed convolution matrix
#    A0_tran = np.zeros([num_theta,num_var])
#    for i in range(num_var): 
#        A0_tran[0,i] = A0[0,i]
#        A0_tran[1::,i] = np.flipud(A0[1::,i])
#    
#    A0ft = np.fft.fft(A0,axis=0)
#    A0_tranft = np.fft.fft(A0_tran,axis=0)    
#    
#    # Enforce Fortran contiguous storage
#    A = np.array(A, order='F', copy=True)
#    A_tran = np.array(A_tran, order='F', copy=True)
#    A0 = np.array(A0, order='F',copy=True)
#    A0_tran = np.array(A0_tran, order='F',copy=True)                                     
#    A0ft = np.array(A0ft, order='F',copy=True)
#    A0_tranft = np.array(A0_tranft, order='F',copy=True)
#
#    # Call Cython function
#    x = cf.circulant_fista(A, A0ft, A0_tran_ft, maxima, b, 
#                           alpha, max_iter, positive=positive)
#        
#    if benchmark: 
#        total_time = time.time() - start_time
#        return x, total_time
#    else:
#        return x    
       
def NPG_ADMM(A,b,beta,lam,eta,rho,max_iters,eps,benchmark=0):
    """
    Inputs:
    A           Matrix
    b           Target vector
    beta        Step size
    max_iters   Max number of iterations
    eps         Convergence threshold
    """
    
    if benchmark:  
        start_time = time.time()
        
    prox_iters = int(max_iters/5)
    
    # Init
    theta = 0
    delta = 10 
    x_old = np.zeros(A.shape[1])  
    x = np.ones(A.shape[1])
       
    # Adaptive step size selection
#    n0 = 0
#    n = 10    
#    xc = 0.50
    
    for i in range(max_iters):
                
        # Iterate
        theta_new = 0.5*(1 + np.sqrt(1 + 4*theta**2))
        x_bar = x + (theta - 1)/theta_new*(x-x_old)
        x_bar[x_bar<0] = 0
        
                
        
        grad = 2*A.T.dot(A.dot(x_bar)-b)
        x_new = prox_ADMM(x_bar - beta*grad,
                          lam,
                          prox_iters,
                          delta,
                          eta,
                          rho)
                          
        delta = norm(x_new - x)/norm(x_new)
        
        # Convergence check
        if delta < eps:
            break      
        
        # Update
        x_old = x
        x = x_new
        

                
        
#        n0 += 1
#        if n0 > n:
#            n0 = 0
#            
#            # Majorization condition
#            diff1 =  (x_new - x)
#            res1 = A.dot(x_new) - b
#            res2 = A.dot(x) - b
#            left = norm(res1)
#            right = norm(res2) + np.dot(diff1,2*A.T.dot(res2)) + 0.5/beta*norm(diff1)**2
#
#            print(left, right)               
#            while( left <= right):
#                beta = beta/xc
#            
#            while( left > right):
#                beta = beta*xc
                
        
    print('---'+str(i)+'---')
    if benchmark: 
        total_time = time.time() - start_time
        return x, total_time
    else:
        return x