# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os
import pickle

#   1. Set up paths: exsitu, data, outputs
# CHESS
# topPath = "/nfs/chess/user/dbanco/VD_sim_processing"
# dataPath = "/nfs/chess/user/seg246/software/MechVD/VD_simulations/c103_polycrystal_sample_2/"
# exsituPath = os.path.join(dataPath,'state_0/simulation/outputs/c103_polycrystal_sample_2_state_0_layer_1_output_data.npz')
# dataFile = os.path.join(dataPath,'state_*/simulation/outputs/c103_polycrystal_sample_2_state_*_layer_1_output_data.npz')

# LOCAL
# topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
topPath = r"E:\VD_sim_processing"
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
params['peak_func'] = "gaussian_rot"
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,'c103_eiger_calibration.yml') #mruby_0401_eiger_calibration
params['roiSize'] = [20,30]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False
params['benchmarkFlag'] = True

# %% 4. Inspect spot tracks on ex-situ and initial scan data
spotInd = 0
scan0 = 0
ttPath = os.path.join(topPath,'outputs_12_13_rot')

scanRange = np.arange(5)
spotsFiles = []
for state in scanRange:
    spotsDir = os.path.join(topPath,f'state_{state}')
    spotsFiles.append(os.path.join(spotsDir,"spots.npz"))

# # %% Analyze track data
# spotInds = [153,268,352,389,585,655,1745]
spotInds = nanInds
resultsData = sf.trackingResultsSim(spotInds,spotsFiles,scanRange,ttPath,params)

# # # Serialize and save the data
# filename = os.path.join(topPath,"results_rot_20.pkl")
# with open(filename, "wb") as file:  # Open file in binary write mode
#     pickle.dump(resultsData, file)

# with open(filename, "rb") as file:  # Open file in binary write mode
#     resultsData = pickle.load(file)

# %% Sort by truthDist
# sumTruthDist = np.nansum(resultsData['truthDist'],1)
# spotIndsSort = np.argsort(sumTruthDist)[::-1]
# sortTruthDist = sumTruthDist[spotIndsSort]

# %% Make track image files

output_path = os.path.join(topPath,'imageFigs_rot_20_30')
dome = 3
num_cols = 5
sf.makeTrackImages(dome,num_cols,output_path,spotInds,spotData,scanRange,dataFile,ttPath,spotsFiles,params)

# %% Make image/Data table slides
# Reorder these so worst tracks are first
ppt_file = os.path.join(output_path, "track_slides.pptx")
sf.makeTrackSlides(ppt_file,output_path,spotInds,resultsData,spotData)

# %% Compute fraction crect detections
# total = len(resultsData['truthDist'].ravel())
# correctDetects = np.sum(resultsData['truthDist'] < 2e-3)
# nanDetects = np.sum(resultsData['truthDist'] == np.nan)

# %% Check continuity of spot truth
spotInds = np.arange(2000,50000)
spotIndArray = sf.gatherSimTruth(spotsFiles,spotInds)
spotIndArray.shape
np.sum(np.isnan(np.sum(spotIndArray,1)))

# nanInds = np.where(np.isnan(np.sum(spotIndArray,1)))[0]
# filename = os.path.join(topPath,"nanInds_50000.pkl")
# with open(filename, "wb") as file:  # Open file in binary write mode
#     pickle.dump(nanInds, file)

