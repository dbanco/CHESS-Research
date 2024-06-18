# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf

# %% Processing Setup

# Output data path
topPath = "C:\\Users\\dpqb1\\Documents\\Data\\ti-2\\"
# Raw HEXD Data
exsituPath = "D:\\CHESS_raw_data\\ti-2-exsitu\\12\\ff\\"
dataPath = "D:\\CHESS_raw_data\\ti-2-tension\\"

# Load in or collect spots data
spotsDir = "spots_11032023"
# Get data from spots data
# sf.collectSpotsData(dataPath, spotsDir)
spotsFile = spotsDir + ".npz"  
spotData = np.load(topPath + spotsFile)
spotInds = np.arange(40,) 


# Full dexela image size and roi size
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = "C:\\Users\\dpqb1\\Documents\\Data\\ti-2\\dex-refined-1.yml"
params['imSize'] = (4888,7300) 
params['roiSize'] = 40
params['gamma'] = [3,5,4,4]

scan1 = 43

sf.spotTracker(dataPath,topPath,exsituPath,spotData,spotInds,params,scan1)
