# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf

# %% Processing Setup

# Output data path
topPath = "/nfs/chess/user/dbanco/"
outPath = "/nfs/chess/user/dbanco/outputs/"

# Raw HEXD Data
exsituPath = "/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-exsitu/12/ff/"
dataPath = "/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/"

# Load in or collect spots data
spotsDir = "spots_11032023"
# Get data from spots data
# sf.collectSpotsData(dataPath, spotsDir)
spotsFile = spotsDir + ".npz"  
spotData = np.load(topPath + spotsFile)
spotInds = np.arange(0,5) 


# Full dexela image size and roi size
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = "/nfs/chess/user/dbanco/dex-refined-1.yml"
params['imSize'] = (4888,7300) 
params['roiSize'] = 40
params['gamma'] = [3,5,3,3]

scan1 = 43

sf.spotTrackerJobs(dataPath,outPath,exsituPath,spotData,spotInds,params,scan1)
