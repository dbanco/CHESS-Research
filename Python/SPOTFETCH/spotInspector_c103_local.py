# -*- coding: utf-8 -*-
"""
Created on Thu May  2 10:07:45 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
import os

# Output data path
# outPath = "C:\\Users\\dpqb1\\Documents\\Data\\c103_2024\\"



# Load in or collect spots data
# topPath = '/nfs/chess/aux/user/wkk32/C103-spots-grains-files'
# spotsDir = 'Sample-{i}/c103-{i}-reconstruction_grains-layer-{j}'
# Get data from spots data
# sf.collectSpotsData(outPath, spotsDir)

topPath = "C:\\Users\\dpqb1\\Documents\\Data\\c103_2024"
spotsDir = "C:\\Users\\dpqb1\\Documents\\Data\\c103_processing"
spotData = np.load(os.path.join(spotsDir,'spots.npz'))

# spotsFile = spotsDir + ".npz"  
# spotData = np.load(outPath + spotsFile)

# Full dexela image size and roi size
params = {}
params['detector'] = 'eiger'
# params['detector'] = 'dexela'
# params['yamlFile'] = "C:\\Users\\dpqb1\\Documents\\Data\\c103_2024\\dexelas_calibrated_ce02_80725_v01.yml"
# params['yamlFile'] = "C:\\Users\\dpqb1\\Documents\\Data\\c103_2024\\dexelas_calibrated_ruby_0504_v01.yml"
params['yamlFile'] = "C:\\Users\\dpqb1\\Documents\\Data\\c103_processing\\eiger16M_monolith_mruby_062224_FINAL.yml"
# params['imSize'] = (4300,7000) 
params['imSize'] = (5000,5000) 
params['roiSize'] = [40,40]   

# interpDirFile = folder_path + "_interp\\" + folder_path+'_interp_frame_'

# %%
# interpDirFile = folder_path + "_interp\\" + folder_path+'_interp_frame_'

# FWHM_etas = np.zeros(numSpots)
# computeTimes = np.zeros(numSpots) 

frame = 100
fnames = "D:\\c103-1-ff-1_0043_EIG16M_CdTe_000383.h5"
# Show HEXD frame with spot wedges
sf.plotSpotWedges(spotData,fnames,frame,params)

# Extract ROIS
# roi_list = sf.loadSpotsAtFrame(spotData,fname1,fname2,frame,params)

# Plot ROIs
sf.plotROIs(roi_list,num_cols = 5)

# %% Plot range of spots
# Notes:
#   ind = 0, frame =4, timeFrms = range(53,68), omeFrms = range(2,7)

# ind = 0
# tth = spotData['tths'][ind]
# eta = spotData['etas'][ind]
# frame = 4
# timeFrms = range(51,61)
# omeFrms = range(-2,9)

# sf.roiAdjacent(yamlFile,ind,tth,eta,imSize,roi_size,frame,omeFrms,timeFrms)

            