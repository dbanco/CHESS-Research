# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 08:53:32 2024

@author: dpqb1
"""
import numpy as np
import spotfetch as sf
import os

topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
dataFile = os.path.join(topPath,r"state_*\simulation\outputs\c103_polycrystal_sample_2_state_*_layer_1_output_data.npz")
dataDir = r"D:\Data\c103"
state = 0
spotsDir = os.path.join(topPath,f'state_{state}','simulation','outputs')
spotsOut = os.path.join(topPath,f'state_{state}')
spotData = np.load(os.path.join(spotsOut,'spots.npz'))

params = {}
params['detector'] = 'eiger_sim'
params['peak_func'] = "gaussian"
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,'c103_eiger_calibration.yml')
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False

grains = [1]
spotInds = sf.findSpots(spotData,grains=grains)

ttPath = os.path.join(topPath,'outputs_1')
fullscanRange = np.arange(0,5)
sf.createTruth(spotInds[6:8],spotData,fullscanRange,ttPath,ttPath,dataFile,params)