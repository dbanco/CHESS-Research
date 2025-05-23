# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 08:53:32 2024

@author: dpqb1
"""
import numpy as np
import spotfetch as sf
import os

# topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
topPath = r"E:\VD_sim_processing"
dataFile1 = os.path.join(topPath,r"state_*\simulation\outputs\c103_polycrystal_sample_2_state_*_layer_1_output_data.npz")
dataFile = os.path.join(topPath,r"state_{t}\simulation\outputs\c103_polycrystal_sample_2_state_{t}_layer_1_output_data.npz")

scanRange = np.arange(5)
spotsFiles = []
for state in scanRange:
    spotsDir = os.path.join(topPath,f'state_{state}')
    spotsFiles.append(os.path.join(spotsDir,"spots.npz"))

params = {}
params['detector'] = 'eiger_sim'
params['peak_func'] = "gaussian_rot"
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,'c103_eiger_calibration.yml') #mruby_0401_eiger_calibration
params['roiSize'] = [30,30]
params['gamma'] = [4,5,9,6,1,1] #[eta,tth,fwhm_eta,fwhm_tth,eta,tth] last two for comparing across omega
params['pool'] = 16
params['parallelFlag'] = False
params['benchmarkFlag'] = True

spotInd = 0

ttPath = os.path.join(topPath,'outputs_1_22_25_rot')
if not os.path.exists(ttPath):
    os.mkdir(ttPath)

dataFileSequence = [] 
for t in scanRange:
    dataFileSequence.append(dataFile.format(t=t)) 

# # Run tracking
# spotInds = np.arange(10957,200000)
# for spotInd in spotInds:
#     initSpotData = np.load(spotsFiles[0])
    # sf.trackSpot(spotInd,initSpotData,dataFileSequence,ttPath,params)

# Process tracking results
print('Processing Results')
spotInds = np.arange(5)
resultsData = sf.trackingResultsSim(spotInds,spotsFiles,scanRange,ttPath,params)
correctDetects = resultsData['truthDist'] < 5
