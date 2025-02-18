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
from scipy.ndimage import center_of_mass

def find_bounding_box(mask):
  """
  Finds the bounding box coordinates of non-zero pixels in a binary mask.

  Args:
    mask: A 3D numpy array representing the binary mask.

  Returns:
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

class Measurememt:
    """Represents a measured blob (candidate spot) in 3D space."""
    
    def __init__(self, blob):
        """
        Initialize a measured spot.
        
        Parameters:
        - com (array): center of mass (tt,eta,ome) of the blob.
        - bound_box (array): boudning box (tth1,tth2,eta1,eta2,ome1,ome2) of blob 
        - intensity (float): Total intensity of the detected spot.
        """
        self.com = center_of_mass(blob)
        self.bound_box = find_bounding_box(blob)
        self.intensity = np.sum(blob)
        
        pass
