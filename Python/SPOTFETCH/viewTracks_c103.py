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
topPath = r"E:\Data\c103_processing"
dataDir = r"E:\Data\c103"

exsituPath = os.path.join(dataDir,"c103-1-ff-1_0004_EIG16M_CdTe_000364.h5")
dataFile = os.path.join(dataDir,"c103-1-ff-1_*_EIG16M_CdTe_{num2:0>6}.h5")
# %% 2. Load in or collect spots data
# Spots for each state file

# sf.collectSpotsData(spotsOut, spotsDir) 
spotData = np.load(os.path.join(topPath,'spots.npz'))

# %% 3. Detector and tracking parameters
params = {}
params['detector'] = 'eiger'
params['peak_func'] = "gaussian_rot"
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,"eiger16M_monolith_mruby_062224_FINAL.yml")
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False
params['benchmarkFlag'] = True

# %% 4. Inspect spot tracks on ex-situ and initial scan data
ttPath = os.path.join(topPath,'outputs_test')
scanRange = np.concatenate((np.array([364,368,372,376,380]), np.arange(383,406), [407]))

# # %% Analyze track data
# resultsData = sf.trackingResultsSim(spotInds,spotsFiles,scanRange,ttPath,params)

# # # Serialize and save the data
# filename = os.path.join(topPath,"results_rot_20.pkl")
# with open(filename, "wb") as file:  # Open file in binary write mode
#     pickle.dump(resultsData, file)

# with open(filename, "rb") as file:  # Open file in binary write mode
#     resultsData = pickle.load(file)


# %% Make track image files
spotInds = [0,1]
output_path = os.path.join(topPath,'imageFigs_c103')
dome = 3
num_cols = 8
sf.makeTrackImages(dome,num_cols,output_path,spotInds,spotData,scanRange,dataFile,ttPath,[],params)

# %% Make image/Data table slides
# ppt_file = os.path.join(output_path, "track_slides_c103.pptx")
# sf.makeTrackSlides(ppt_file,output_path,spotInds,resultsData,spotData)

# %% Compute fraction crect detections
# total = len(resultsData['truthDist'].ravel())
# correctDetects = np.sum(resultsData['truthDist'] < 2e-3)
# nanDetects = np.sum(resultsData['truthDist'] == np.nan)
