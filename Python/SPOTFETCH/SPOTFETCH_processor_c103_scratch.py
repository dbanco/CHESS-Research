# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
import os

# %% Processing Setup

# Output data path
topPath = "/mnt/scratch/dbanco/c103_processing/Sample-1/layer-1"

# Raw HEXD Data
exsituPath = "/nfs/chess/raw/2024-2/id3a/miller-3528-c/c103-1-ff-1"
dataPath = "/mnt/scratch/....TBD"

# Load in or collect spots data
spotsPath = os.path.join(topPath,'spots','spots.npz')

# Get data from spots data  
spotData = np.load(spotsPath)
spotInds = np.arange(0,25) 

# Full dexela image size and roi size
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'eiger'
params['yamlFile'] = '/mnt/scratch/dbanco/mruby_eiger_calibration_single_grain_v01.yml'
params['imSize'] = (5000,5000) 
params['roiSize'] = 40
params['gamma'] = [8,8,4,4]

scan1 = 1

sf.spotTracker(dataPath,topPath,exsituPath,spotData,spotInds,params,scan1)
