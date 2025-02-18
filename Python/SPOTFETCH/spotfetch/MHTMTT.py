# -*- coding: utf-8 -*-
"""
spotfetch.MHTMTT

Description:
Module for tracking x-ray diffraction spots in 3D using mulitple hypotheses test
multitarget tracking
Created on Tue Feb 18 09:59:02 2025

@author: dpqb1
"""

import numpy as np
from scipy.spatial.distance import mahalanobis
from scipy.optimize import linear_sum_assignment
from scipy.ndimage import gaussian_filter
from scipy.ndimage import center_of_mass
from skimage.feature import hessian_matrix
from scipy.ndimage import label

class Measurement:
    """Represents a measurement (candidate spot) in 3D space."""
    
    def __init__(self, x):
        """
        Initialize a measured spot.
        
        Parameters:
        -----------
            - x (array): masked pixel intesity data
        
        Properties:
        - com (array): center of mass (tt,eta,ome) of the blob.
        - bound_box (array): boudning box (tth1,tth2,eta1,eta2,ome1,ome2) of blob 
        - intensity (float): Total intensity of the detected spot.
        """
        self.com = center_of_mass(x)
        self.bound_box = find_bounding_box(x)
        self.intensity = np.sum(x)
        
        pass

def find_bounding_box(mask):
  """
  Finds the bounding box coordinates of non-zero pixels in a binary mask.

  Parameters:
  -----------
    mask: A 3D numpy array representing the binary mask.

  Returns:
  --------
    A tuple (tth_min, tth_max, eta_min, eta_max, ome_min, ome_max) representing 
    the bounding box coordinates or None if no non-zero pixels are found.
  """
  rows = np.any(mask, axis=1)
  cols = np.any(mask, axis=0)
  tubs = np.any(mask, axis=2)
  
  if not np.any(rows) or not np.any(cols) or not np.any(tubs):
    return None  # Return None if the mask is empty

  tth_min, tth_max = np.where(rows)[0][[0, -1]]
  eta_min, eta_max = np.where(cols)[0][[0, -1]]
  ome_min, ome_max = np.where(tubs)[0][[0, -1]]
  
  return tth_min, tth_max, eta_min, eta_max, ome_min, ome_max

def detectBlobHDoG(x,params):
    """
    Detects blobs using the Hessian-based Difference of Gaussians (DoG) method.

    Parameters:
    -----------
    x : ndarray
        Input 3D data.

    Returns:
    --------
    tuple
        - blobs : ndarray
            Labeled blob regions.
        - num_blobs : int
            Number of blobs detected.
        - hess_mat : list of ndarray
            Hessian matrices of the DoG result.
    """
    # 1. Compute normalized DoG
    dog_norm = DoG(x, sigma=params['sigmas'], dsigma=params['dsigma'])

    # 2. Pre-segmentation
    hess_mat = hessian_matrix(dog_norm,use_gaussian_derivatives=params['gaussian_derivatives'])
    D1 = np.zeros(hess_mat[0].shape)
    D2 = np.zeros(hess_mat[0].shape)
    D3 = np.zeros(hess_mat[0].shape)
    for i1 in range(hess_mat[0].shape[0]):
        for i2 in range(hess_mat[0].shape[1]):
            for i3 in range(hess_mat[0].shape[2]):
                h_mat = np.array([
                    [hess_mat[0][i1, i2, i3], hess_mat[1][i1, i2, i3], hess_mat[2][i1, i2, i3]],
                    [hess_mat[1][i1, i2, i3], hess_mat[3][i1, i2, i3], hess_mat[4][i1, i2, i3]],
                    [hess_mat[2][i1, i2, i3], hess_mat[4][i1, i2, i3], hess_mat[5][i1, i2, i3]]
                ])
                D1[i1,i2,i3] = h_mat[0,0]
                D2[i1,i2,i3] = np.linalg.det(h_mat[:2,:2])
                D3[i1,i2,i3] = np.linalg.det(h_mat)

    posDefIndicator = (D1 > 0) & (D2 > 0) & (D3 > 0)
    blobs, num_blobs = label(posDefIndicator)
    return blobs, num_blobs, hess_mat

def DoG(x, sigma, dsigma, gamma=2):
    """
    Computes the Difference of Gaussians (DoG) approximation.

    Parameters:
    -----------
    x : ndarray
        Input data.
    sigma : float
        Base Gaussian sigma.
    dsigma : float
        Increment for sigma.
    gamma : float, optional
        Scaling factor, default is 2.

    Returns:
    --------
    ndarray
        Normalized DoG result.
    """
    g1 = gaussian_filter(x, sigma=sigma)
    g2 = gaussian_filter(x, sigma=sigma + dsigma)
    return (g2 - g1) / (np.mean(sigma) * np.mean(dsigma))

def blobsToMeasurements(x,blobs):
    """
    Converts labeled blobs into a list of Measurement objects.

    Parameters:
    -----------
    blobs : ndarray
        Labeled blob regions.

    Returns:
    --------
    measurements:
        List of measurement objects containing candidate spots.
    """
    measurements = []
    unique_labels = np.unique(blobs)
    
    # Exclude background (label 0)
    unique_labels = unique_labels[unique_labels != 0]

    for blob_label in unique_labels:
        mask = (blobs == blob_label)  # Extract binary mask for the current blob
        x_masked = x * mask
        measurement = Measurement(x_masked)
        measurements.append(measurement)

    return measurements