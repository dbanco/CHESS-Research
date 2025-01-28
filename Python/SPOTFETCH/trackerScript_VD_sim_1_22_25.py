# -*- coding: utf-8 -*-
"""
Runs 2D Gaussian tracking on VD simulated data

Created on Tue Nov 12 08:53:32 2024
Author: Daniel Banco
Email: dpqb10@gmail.com
"""
import numpy as np
import spotfetch as sf
import os

# topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
topPath = r"E:\VD_sim_processing"
dataFile1 = os.path.join(topPath,r"state_*\simulation\outputs\c103_polycrystal_sample_2_state_*_layer_1_output_data.npz")
dataFile = os.path.join(topPath,r"state_{scan}\simulation\outputs\c103_polycrystal_sample_2_state_{scan}_layer_1_output_data.npz")

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
params['gamma'] = [4,5,9,6,2,2] #[eta,tth,fwhm_eta,fwhm_tth,eta,tth] last two for comparing across omega
params['pool'] = 16
params['parallelFlag'] = False
params['benchmarkFlag'] = True

prmBlob = {}
prmBlob['detector'] = 'eiger_sim'
prmBlob['imSize'] = (5000,5000)
prmBlob['yamlFile'] = os.path.join(topPath,'c103_eiger_calibration.yml') #mruby_0401_eiger_calibration
prmBlob['roiSize'] = [31,31,11]
prmBlob['gamma'] = [4,5,9,6,2,2]
prmBlob['pool'] = 16
prmBlob['parallelFlag'] = False
prmBlob['benchmarkFlag'] = True

spotInd = 0
suffix = '_1_22_25_rot_chess'
ttPath = os.path.join(topPath,'outputs'+suffix)

if not os.path.exists(ttPath):
    os.mkdir(ttPath)

dataFileSequence = [] 
for scan in scanRange:
    dataFileSequence.append(dataFile.format(scan=scan)) 
    initSpotData = np.load(spotsFiles[0])

spotInds = np.arange(8000)

# # Run tracking
# for spotInd in spotInds:
#     sf.trackSpot(spotInd,initSpotData,dataFileSequence,ttPath,params)



# Process tracking results
print('Processing Results')
resultsData = sf.trackingResultsSim(spotInds,spotsFiles,scanRange,ttPath,params)
correctDetects = resultsData['truthDist'] < 5

# Make track image files
# output_path = os.path.join(topPath,'imageFigs'+suffix)
# dome = 4
# num_cols = 5
# sf.makeTrackImages(dome,num_cols,output_path,spotInds,initSpotData,scanRange,dataFile,ttPath,spotsFiles,params)

