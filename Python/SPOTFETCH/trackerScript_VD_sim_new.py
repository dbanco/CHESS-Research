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
params['roiSize'] = [20,30]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False
params['benchmarkFlag'] = True

spotInd = 0

ttPath = os.path.join(topPath,'outputs_12_13_rot')
if not os.path.exists(ttPath):
    os.mkdir(ttPath)

dataFileSequence = [] 
for t in scanRange:
    dataFileSequence.append(dataFile.format(t=t)) 

spotInds = np.arange(100)
for spotInd in nanInds:
    initSpotData = np.load(spotsFiles[0])
    sf.trackSpot(spotInd,initSpotData,dataFileSequence,ttPath,params)

# Update this function so it's consistent with trackSpot maybe?
# print('Processing Results')
# resultsData = sf.trackingResultsSim(spotInds,spotsFiles,scanRange,ttPath,params)
# correctDetects = resultsData['truthDist'] < 5
