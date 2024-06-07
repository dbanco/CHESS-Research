# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 09:22:38 2024

@author: dpqb1
"""

import numpy as np
# import matplotlib.pyplot as plt
# import os
import spotfetch as sf
# from hexrd.fitting import fitpeak
# from hexrd.fitting import peakfunctions as pkfuncs
# import time
# import h5py

# %% First load indexed spot data
folder_path = "spots_11032023"  
spot_data = np.load(folder_path+'.npz')    

# %% Dataset Parameters
fDir = "D:\\CHESS_raw_data\\ti-2-exsitu\\"
fname1 = fDir + '12\\ff\\ff1_000098.h5'
fname2 = fDir + '12\\ff\\ff2_000098.h5'
yamlFile = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\dex-refined-1.yml"
imSize = [4888,7300]
center = [imSize[0]/2,imSize[1]/2]
roi_size = 40

# %%
interpDirFile = folder_path + "_interp\\" + folder_path+'_interp_frame_'
numSpots = 46181
# FWHM_etas = np.zeros(numSpots)
# computeTimes = np.zeros(numSpots) 

dirNum = 1
fileNum = 105 + dirNum
frame = 39
dNum = str(dirNum)
fNum = str(fileNum)
fDir = "D:\\CHESS_raw_data\\ti-2-tension\\"
fname1 = fDir + dNum + '\\ff\\ff1_000' + fNum + '.h5'
fname2 = fDir + dNum + '\\ff\\ff2_000' + fNum + '.h5'

# Extract ROIS
roi_list = sf.loadSpotsAtFrame(spot_data,fname1,fname2,yamlFile,frame,imSize,roi_size)

# Plot ROIs
sf.plotROIs(roi_list,num_cols = 5)

# Plot ROIs on full image
# sf.plotSpotWedges(spot_data,fname1,fname2,yamlFile,frame,center,roi_size)
    
# sf.plotSpotRectangles(spot_data,fname1,fname2,yamlFile,frame,center,roi_size)