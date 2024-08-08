# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 10:19:33 2024

@author: dpqb1
"""
import spotfetch as sf
import numpy as np
import os

##### Paths
topPath = "/nfs/chess/user/dbanco/ti-2_processing"
dataPath = "/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/"
trackPath = os.path.join(topPath,'outputs')

# Load in or collect spots data
spotsDir = "spots_11032023"
# Get data from spots data
# sf.collectSpotsData(dataPath, spotsDir)
spotsFile = spotsDir + ".npz"  
spotData = np.load(os.path.join(topPath,'spots',spotsFile))

# Full dexela image size and roi size
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = "/nfs/chess/user/dbanco/ti-2_processing/dex-refined-1.yml"
params['imSize'] = (4888,7300) 
params['roiSize'] = 40

dome = 4
scanRange = np.arange(43,62)
spotInds = [8]

sf.roiTrackVisual(spotInds,spotData,dome,scanRange,trackPath,dataPath,params)

