# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os


#   1. Set up paths: exsitu, data, outputs
topPath = "/nfs/chess/user/dbanco/VD_sim_processing"
dataPath = "/nfs/chess/user/seg246/software/MechVD/VD_simulations/c103_polycrystal_sample_2/"

# exsituPath = os.path.join(dataPath,f'c103-1-ungripped-1_{num1:0>4}_EIG16M_CdTe_000353.h5')
exsituPath = os.path.join(dataPath,'state_0/simulation/outputs/c103_polycrystal_sample_2_state_0_layer_1_output_data.npz')
dataFile = os.path.join(dataPath,'state_*/simulation/outputs/c103_polycrystal_sample_2_state_*_layer_1_output_data.npz')


# %% 2. Load in or collect spots data
# Spots for each state file
state = 0
spotsDir = os.path.join(dataPath,f'state_{state}/simulation/outputs')
spotsOut = os.path.join(topPath,f'state_{state}')
    
# sf.collectSpotsData(spotsOut, spotsDir) 
spotData = np.load(os.path.join(spotsOut,'spots.npz'))

# %% 3. Detector and tracking parameters
params = {}
params['detector'] = 'eiger_sim'
params['peak_func'] = 'Gaussian'
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(dataPath,'c103_eiger_calibration.yml')
# params['detector'] = 'dexela'
# params['imSize'] = (4888,7300) 
# params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16

# %% 4. Inspect spot tracks on ex-situ and initial scan data
grains = [44,158]
# tth110 = 0.06562438
# tth200 = 0.09267698
# tth211 = 0.1136209
# tth321 = 0.1736603
# tths = [tth110, tth211]
# hklNames = ['110', '211']
# spotInds = sf.findSpots(spotData,grains=grains, tth=tths, dtth=0.012)
spotInds = sf.findSpots(spotData,grains=grains)
# spotInds = [113,205,413,801]


sf.plotSpotWedges(spotData,exsituPath,100,params)

scanRange = np.arange(1,6)
trackPath = os.path.join(topPath,'outputs')

# # %% 6. Begin Processing
initTracksPath = os.path.join(topPath,'outputs')
# sf.initExsituTracks(initTracksPath,exsituPath,spotData, spotInds, params, 0)

# advance = False
# num1 = 1
# for num2 in scanRange:
#     sf.spotTracker(dataFile,topPath,spotData,spotInds,params,num1,num2,advance)

# spotInds = [113,205,413,801]
dome = 2
scanRange = np.arange(0,6)
sf.roiTrackVisual(spotInds,spotData,dome,scanRange,trackPath,dataFile,params)