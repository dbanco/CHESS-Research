# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
# import matplotlib.pyplot as plt
# import pickle
# import os

# %% Processing Setup

# Output data path
topPath = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains"
output_path = topPath + "\\outputs"
output_file = "trackData.pkl"
outFile = output_path + output_file

# Initial indexing data
folder_path = "spots_11032023"  
spotData = np.load(folder_path+'.npz')
# spotInds = np.arange(0,4)
# frame = 4
# spotInds = np.where(spot_data['ome_idxs'] == frame)[0]
grain = 0
spotInds = np.where(spotData['id_nums'] == grain)[0]

# HEXD Data
dataPath = "D:\\CHESS_raw_data\\ti-2-tension\\"

# Full dexela image size and roi size
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\dex-refined-1.yml"
params['imSize'] = (4888,7300) 
params['roiSize'] = 40

scan1 = 43

sf.spotTracker(dataPath,outFile,spotData,spotInds,params,scan1)
