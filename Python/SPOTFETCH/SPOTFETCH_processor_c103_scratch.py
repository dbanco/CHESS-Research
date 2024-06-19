# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
import os
import shutil

# %% Processing Setup

# Output data path
topPath = "/nfs/chess/user/dbanco/c103_processing/Sample-1/layer-1"
# topPath = "/mnt/scratch/dbanco/c103_processing/Sample-1/layer-1"

# Raw HEXD Data
exsituPath = "/nfs/chess/raw/2024-2/id3a/miller-3528-c/c103-1-ff-1/1/ff"
dataPath = "/mnt/scratch/....TBD"

# Load in or collect spots data
spotsPath = os.path.join(topPath,'spots','spots.npz')

# Get data from spots data  
spotData = np.load(spotsPath)
spotInds = np.arange(0,25) 

# Full dexela image size and roi size
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
# Exsitu Yaml
params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
# params['yamlFile'] = '/mnt/scratch/dbanco/mruby_eiger_calibration_single_grain_v01.yml'
params['imSize'] = (4888,7300) 
# params['imSize'] = (5000,5000) 
params['roiSize'] = 40
params['gamma'] = [8,8,4,4]

scan1 = 1
scan0 = scan1 - 1 
# Run on CHESS machine
initTracksPath = os.path.join(topPath,'initTracks')
sf.initExsituTracks(initTracksPath,exsituPath, spotData, spotInds, params, scan0)

# # Run on Wrangler
# Copy track files from init to outputs
# for k in spotInds:
#     src = os.path.join(initTracksPath,'trackData_{k}.pkl')
#     des = os.path.join(topPath,'outputs')
#     shutil.copyfile(src,des)

# sf.spotTracker(dataPath,topPath,exsituPath,spotData,spotInds,params,scan1)