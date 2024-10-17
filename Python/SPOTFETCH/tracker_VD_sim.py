# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os

dataFile = "C:\\Users\\dpqb1\\Documents\\Data\\c103_polycrystal_sample_2_state_0_layer_1_output_data.npz"
simData = np.load(dataFile)
shp = simData['shape']
frm = 0
simImg = np.zeros((shp[0],shp[1]))

rowD = simData[f'{frm}_row']
colD = simData[f'{frm}_col']
datD = simData[f'{frm}_data']

for i in range(len(rowD)):
    simImg[rowD[i],colD[i]] = datD[i]
    
# plt.imshow(simImg)
# plt.colorbar()
# plt.clim((0,10))
    
#   1. Set up paths: exsitu, data, outputs
topPath = "/nfs/chess/user/dbanco/VD_sim_processing"
dataPath = "/nfs/chess/user/seg246/software/MechVD/VD_simulations/c103_polycrystal_sample_2/"

# exsituPath = os.path.join(dataPath,f'c103-1-ungripped-1_{num1:0>4}_EIG16M_CdTe_000353.h5')
exsituPath = os.path.join(dataPath,"state0/simulation/outputs/c103_polycrystal_sample_2_state_0_layer_1_output_data.npz")
dataFile = os.path.join(dataPath,"state*/simulation/outputs/c103_polycrystal_sample_2_state_*_layer_1_output_data.npz")


# %% 2. Load in or collect spots data
# Spots for each state file
spotsDir = os.path.join(dataPath,"state*/simulation/outputs/eiger")
spotsOut = os.path.join(topPath,'state*')

sf.collectSpotsData(spotsOut, spotsDir) 
spotData = np.load(os.path.join(spotsOut,'spots.npz'))

# %% 3. Detector and tracking parameters
params = {}
params['detector'] = 'eiger'
params['peak_func'] = 'Gaussian'
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,'c103_eiger_calibration.yml')
# params['detector'] = 'dexela'
# params['imSize'] = (4888,7300) 
# params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]0
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

scanRange = np.arange(1,6)
trackPath = os.path.join(topPath,'outputs')

# # %% 6. Begin Processing
initTracksPath = os.path.join(topPath,'outputs')
# sf.initExsituTracks(initTracksPath,exsituPath,spotData, spotInds, params, 364)

advance = False
num1 = 1
for num2 in scanRange:
    sf.spotTracker(dataFile,topPath,spotData,spotInds,params,num1,num2,advance)

# spotInds = [113,205,413,801]
# dome = 3
# sf.roiTrackVisual(spotInds,spotData,dome,scanRange,trackPath,dataFile,params)