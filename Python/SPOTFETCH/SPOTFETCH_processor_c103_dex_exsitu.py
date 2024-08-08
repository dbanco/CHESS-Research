# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf

# %% Processing Setup

# Output data path
outPath = "C:\\Users\\dpqb1\\Documents\\Data\\c103_2024\\"

# Track file
trackFile = "trackData.pkl"
trackPath = outPath + "outputs\\" + trackFile

# Raw HEXD Data
dataPath = "C:\\Users\\dpqb1\\Documents\\Data\\c103_2024\\c103-1-ff-1\\"

# Load in or collect spots data
spotsDir = "Sample-4\\c103-4-reconstruction_grains-layer-1"
# Get data from spots data
# sf.collectSpotsData(outPath, spotsDir)
spotsFile = spotsDir + ".npz"  
spotData = np.load(outPath + spotsFile)
# spotInds = np.arange(0,4)
# frame = 4
# spotInds = np.where(spot_data['ome_idxs'] == frame)[0]
grain = 0
spotInds = np.where(spotData['id_nums'] == grain)[0]

# Full dexela image size and roi size
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = "C:\\Users\\dpqb1\\Documents\\Data\\c103_2024\\dexelas_calibrated_ce02_80725_v01.yml"
params['imSize'] = (4888,7300) 
params['roiSize'] = 40

scan1 = 43

sf.spotTracker(dataPath,trackFile,spotData,spotInds,params,scan1)
