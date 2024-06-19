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
topPath_nfs = "/nfs/chess/user/dbanco/c103_processing/Sample-1/layer-1"
topPath_mnt = "/mnt/scratch/dbanco/c103_processing/Sample-1/layer-1"

# Raw HEXD Data
exsituPath = "/nfs/chess/raw/2024-2/id3a/miller-3528-c/c103-1-ff-1/1/ff"
dataPath = "/mnt/scratch/....TBD"
initTracksPath_nfs = os.path.join(topPath_nfs,'initTracks')
initTracksPath_mnt = os.path.join(topPath_mnt,'initTracks')

# Load in or collect spots data
spotsPath = os.path.join(topPath_nfs,'spots','spots.npz') 
spotData = np.load(spotsPath)
K = len(spotData['etas'])
spotInds = np.arange(0,K) 

# Full dexela image size and roi size
params = {}

params['detector'] = 'dexela'
params['imSize'] = (4888,7300) 
params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'

# params['detector'] = 'eiger'
# params['imSize'] = (5000,5000) 
# params['yamlFile'] = '/mnt/scratch/dbanco/mruby_eiger_calibration_single_grain_v01.yml'

params['roiSize'] = 40
params['gamma'] = [8,8,4,4]

scan1 = 1
scan0 = scan1 - 1 
# Run on CHESS machine
sf.initExsituTracks(initTracksPath_mnt,exsituPath, spotData, spotInds, params, scan0)

# # Run on Wrangler
# # Copy track files from init to outputs
# for k in spotInds:
#     if np.mod(k,1000) == 0:
#         print(k)
#     src = os.path.join(initTracksPath_nfs,f'trackData_{k}.pkl')
#     des = os.path.join(initTracksPath_mnt,f'trackData_{k}.pkl')
#     # des = os.path.join(topPath_mnt,'outputs',f'trackData_{k}.pkl')
#     shutil.copyfile(src,des)

# sf.spotTracker(dataPath,topPath,exsituPath,spotData,spotInds,params,scan1)