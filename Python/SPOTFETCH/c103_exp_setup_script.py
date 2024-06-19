# -*- coding: utf-8 -*-
"""
Created on Thu May  2 10:07:45 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
import os
import matplotlib.pyplot as plt

# Output data path
# Sets up directories and loads in the spot files
for i in range(1,3):
    for j in range(1,3):
        outPath = f"/mnt/scratch/dbanco/c103_processing/Sample-{i}"
        outPath2 = os.path.join(outPath,f'layer-{j}')
        outPath3 = os.path.join(outPath2,'spots')
        outPath4 = os.path.join(outPath2,'jobs')
        outPath5 = os.path.join(outPath2,'outputs')
        outPath6 = os.path.join(outPath2,'initTracks')
        for path_i in [outPath,outPath2,outPath3,outPath4,outPath5,outPath6]:
            if not os.path.exists(path_i):
                os.mkdir(path_i)
                
        # Load in or collect spots data
        topPath = '/nfs/chess/aux/user/wkk32/C103-spots-grains-files'
        spotsDir = f'Sample-{i}/c103-{i}-reconstruction_grains-layer-{j}'
        
        # Get data from spots data
        # sf.collectSpotsData(outPath3, os.path.join(topPath,spotsDir))
      
# %%
fname = '/mnt/scratch/pagan-dwellfatigue-series-cycle-2024-2/ti6242-sample1-2_0011_EIG16M_CdTe_000001.h5'
spotsPath = f"/mnt/scratch/dbanco/c103_processing/Sample-1/layer-1/spots"

params = {}
params['detector'] = 'eiger'
params['yamlFile'] = '/mnt/scratch/dbanco/mruby_eiger_calibration_single_grain_v01.yml'
params['imSize'] = (5000,5000) 
params['roiSize'] = 40  

# Show HEXD frame with spot wedges
frame = 10
spotData = np.load(os.path.join(spotsPath,'spots.npz'))
sf.plotSpotWedges(spotData,fname,fname,frame,params)

# Extract ROIS
# roi_list = sf.loadSpotsAtFrame(spotData,fname1,fname2,frame,params)

# Plot ROIs
# sf.plotROIs(roi_list,num_cols = 5)

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

            