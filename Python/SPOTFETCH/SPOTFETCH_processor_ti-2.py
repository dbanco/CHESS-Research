# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
import os
import numpy as np
import spotfetch as sf

# %% Processing Setup

# Output data path
topPath = "C:\\Users\\dpqb1\\Documents\\Data\\ti-2\\"
# Raw HEXD Data
exsituPath = "D:\\CHESS_raw_data\\ti-2-exsitu\\12\\ff\\"
dataPath = "D:\\CHESS_raw_data\\ti-2-tension\\"

# Load in or collect spots data 
spotData = np.load(os.path.join(topPath,'spots','spots.npz'))
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
