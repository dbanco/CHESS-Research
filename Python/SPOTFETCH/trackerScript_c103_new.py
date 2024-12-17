# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 08:53:32 2024

@author: dpqb1
"""
import numpy as np
import spotfetch as sf
import os

topPath = r"C:\Users\dpqb1\Documents\Data\c103_processing"
dataDir = r"E:\Data\c103"

dataFile = os.path.join(dataDir,"c103-1-ff-1_*_EIG16M_CdTe_{num2:0>6}.h5")
scanRange = np.concatenate(( np.array([364,368,372,376,380]), 
                             np.arange(383,406), [407] ))

# spotsDir = r"C:\Users\dpqb1\Documents\Data\c103_2024\Sample-1\c103-1-reconstruction_grains-layer-1"
spotsDir = r"E:\Data\c103\C103_1_unloaded_gripped_layer2of4\layer_1\3A_grains"
# sf.collectSpotsData(topPath,spotsDir)

spotsFile = os.path.join(topPath,"spots.npz")

spotData = np.load(spotsFile)

params = {}
params['detector'] = 'eiger'
params['peak_func'] = "gaussian_rot"
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,"eiger16M_monolith_mruby_062224_FINAL.yml")
params['roiSize'] = [30,30]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False
params['benchmarkFlag'] = True

grains = [44,158]
spotInds = sf.findSpots(spotData,grains=grains)
# spotInds = np.arange(113)
spotInd = 0

ttPath = os.path.join(topPath,'outputs_12_13')
if not os.path.exists(ttPath):
    os.mkdir(ttPath)

dataFileSequence = sf.getDataFileSequence(dataFile,scanRange)   
sf.trackSpot(spotInd,spotData,dataFileSequence,ttPath,params)

# %%
spotInds = [0,1]
output_path = os.path.join(topPath,'imageFigs_c103')
dome = 3
num_cols = 8
sf.makeTrackImages(dome,num_cols,output_path,spotInds,spotData,scanRange,dataFile,ttPath,[],params)

