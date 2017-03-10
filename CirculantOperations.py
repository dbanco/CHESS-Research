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
from functools import partial

def convolve_2D(aft,bft):
    c = np.fft.ifft2(aft*bft).real
    return c

def convolve_2D_b(a_bft):
    c = np.fft.ifft2(a_bft[0]*a_bft[1]).real
    return c


def convolve_2D_svd(svd_x_list,Rft):
    u   = svd_x_list[0]
    s   = svd_x_list[1]
    v   = svd_x_list[2]
    x   = svd_x_list[3]
    c = np.fft.ifft2(s*np.outer(u,v)*np.fft.fft2(x)).real
    return c

def convolve_2D_svd_b(svd_x_list):
    u   = svd_x_list[0]
    s   = svd_x_list[1]
    v   = svd_x_list[2]
    x   = svd_x_list[3]
    c = np.fft.ifft2(s*np.outer(u,v)*np.fft.fft2(x)).real
    return c

def add_x_to_list(A0ft_list, x):
    k = 0
    A0ft_x_list = A0ft_list[:]
    for t in enumerate(x.shape[2]):
        for r in enumerate(x.shape[3]):
            A0ft_x_list[k].append(x[:,:,t,r])
            k += 1
            return A0ft_x_list

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

"Circulant matrix-vector product subroutine"
def AtR_ft_2D_Parallel(A0ft_list, R):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        R       Residual vector
    """
    
    R_ft = np.fft.fft2(R)

    pool = Pool()
        
    partial_convolve = partial(convolve_2D,bft=R_ft)

    AtR = np.asarray(pool.map(partial_convolve,A0ft_list))  

    print(AtR.shape)
    
    return AtR
    

"Circulant matrix-vector product subroutine"
def Ax_ft_2D_Parallel(A0ft_list, x):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        x       Coefficient vector
    """
    
    A0ft_xft_list = add_x_to_list(A0ft_list,x)
    
    pool = Pool() 
    Ax = np.asarray(pool.map(convolve_2D_b,A0ft_xft_list))
    
    
    return np.sum(Ax,0)
  
    "Circulant matrix-vector product subroutine"
def AtR_ft_2D_Parallel_SVD(A0ft_list, R):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        R       Residual vector
    """
    
    R_ft = np.fft.fft2(R)

    pool = Pool()
        
    partial_convolve = partial(convolve_2D_svd,Rft=R_ft)

    AtR = np.asarray(pool.map(partial_convolve,A0ft_list))  

    print(AtR.shape)
    
    return AtR
    

"Circulant matrix-vector product subroutine"
def Ax_ft_2D_Parallel_SVD(A0ft_list_svd, x):  
    """
    Inputs:
        A0ft    Each column is the first column of a circulant matrix Ai
        x       Coefficient vector
    """
    
    svd_x_list = add_x_to_list(A0ft_list_svd,x)
    
    pool = Pool() 
    Ax = np.asarray(pool.map(convolve_2D_svd_b,svd_x_list))
    
    
    return np.sum(Ax,0)