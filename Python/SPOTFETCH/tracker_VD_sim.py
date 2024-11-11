# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os
import matplotlib.pyplot as plt


#   1. Set up paths: exsitu, data, outputs
# CHESS
# topPath = "/nfs/chess/user/dbanco/VD_sim_processing"
# dataPath = "/nfs/chess/user/seg246/software/MechVD/VD_simulations/c103_polycrystal_sample_2/"
# exsituPath = os.path.join(dataPath,'state_0/simulation/outputs/c103_polycrystal_sample_2_state_0_layer_1_output_data.npz')
# dataFile = os.path.join(dataPath,'state_*/simulation/outputs/c103_polycrystal_sample_2_state_*_layer_1_output_data.npz')

# LOCAL
topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
exsituPath = os.path.join(topPath,r"state_0\simulation\outputs\c103_polycrystal_sample_2_state_0_layer_1_output_data.npz")
dataFile = os.path.join(topPath,r"state_*\simulation\outputs\c103_polycrystal_sample_2_state_*_layer_1_output_data.npz")

# %% 2. Load in or collect spots data
# Spots for each state file
state = 0
# CHESS
# spotsDir = os.path.join(dataPath,f'state_{state}','simulation','outputs')

spotsDir = os.path.join(topPath,f'state_{state}','simulation','outputs')
spotsOut = os.path.join(topPath,f'state_{state}')
    
# sf.collectSpotsData(spotsOut, spotsDir) 
spotData = np.load(os.path.join(spotsOut,'spots.npz'))

# %% 3. Detector and tracking parameters
params = {}
params['detector'] = 'eiger_sim'
params['peak_func'] = "gaussian"
params['imSize'] = (5000,5000)
# params['yamlFile'] = os.path.join(dataPath,'c103_eiger_calibration.yml')
params['yamlFile'] = os.path.join(topPath,'c103_eiger_calibration.yml')
# params['detector'] = 'dexela'
# params['imSize'] = (4888,7300) 
# params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False

# %% 4. Inspect spot tracks on ex-situ and initial scan data
grains = [1]
# tth110 = 0.06562438
# tth200 = 0.09267698
# tth211 = 0.1136209
# tth321 = 0.1736603
# tths = [tth110, tth211]
# hklNames = ['110', '211']
# spotInds = sf.findSpots(spotData,grains=grains, tth=tths, dtth=0.012)
spotInds = sf.findSpots(spotData,grains=grains)
# spotInds = sf.findSpots(spotData)
# spotInds = [113,205,413,801]

# spotData = np.load(os.path.join(spotsOut,'spots.npz'))
# sf.plotSpotWedges(spotData,exsituPath,0,params)

# # testInd = [9228]
# testInd = [9228,20987]
testInd = spotInds

trackPath = os.path.join(topPath,'outputs_1')
os.mkdir(trackPath)


# %% Inspect spots
    
# roi_list = sf.loadSpotsAtFrame(spotData,exsituPath,2,params)
# sf.plotROIs(roi_list,5)


# %% 6.Tracking
sf.initExsituTracks(trackPath,exsituPath,spotData, testInd, params, 0)
advance = False
num1 = 1
scanRange = np.arange(1,5)
for scan in scanRange:
    sf.spotTracker(dataFile,trackPath,spotData,testInd,params,num1,scan,advance)

# 7. True spot locations for all states
# spotFiles = []
# fullscanRange = np.arange(0,5)
# for state in fullscanRange:
#     spotsOut = os.path.join(topPath,f'state_{state}')
#     spotFiles.append(os.path.join(spotsOut,'spots.npz'))
    
# # # 8. Tracking results
# trackPath = os.path.join(topPath,'outputs_1')
# foundInd, numTracks0, changeDist, truthDist = sf.trackingResults(testInd,spotFiles,fullscanRange,trackPath,dataFile,params)


# plt.hist(truthDist)
# plt.hist(changeDist[:,1:5])


# 9. Show tracks
# dome = 3
# sf.roiTrackVisual(testInd[6:10],spotData,dome,fullscanRange,trackPath,dataFile,params)