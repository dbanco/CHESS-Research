# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
import os
import time

# %% Processing Setup

# Output data path
topPath = "/nfs/chess/user/dbanco/ti-2_processing"

# Raw HEXD Data
exsituPath = "/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-exsitu/12/ff/"
dataPath = "/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/"

# Load in or collect spots data
spotsDir = "spots_11032023"
# Get data from spots data
# sf.collectSpotsData(dataPath, spotsDir)
spotsFile = spotsDir + ".npz"  
spotData = np.load(os.path.join(topPath,'spots',spotsFile))
# spotInds = sf.findSpots(spotData,5,np.pi/2,0.1)
# spotInds = np.arange(120,1000)
spotInds = np.arange(8)

# Full dexela image size and roi size
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = os.path.join(topPath,'dex-refined-1.yml')
params['imSize'] = (4888,7300) 
params['roiSize'] = 40
params['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]

scan1 = 43
initTracksPath = os.path.join(topPath,'outputs')
sf.initExsituTracks(initTracksPath,exsituPath, spotData, spotInds, params, scan1-1)
sf.spotTracker(dataPath,topPath,exsituPath,spotData,spotInds,params,scan1)
