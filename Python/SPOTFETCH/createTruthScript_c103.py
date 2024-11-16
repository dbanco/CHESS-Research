# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 08:53:32 2024

@author: dpqb1
"""
import numpy as np
import spotfetch as sf
import os
import tkinter as tk

topPath = r"C:\Users\dpqb1\Documents\Data\c103_processing"
dataDir = r"E:\Data\c103"

exsituPath = os.path.join(dataDir,"c103-1-ff-1_0004_EIG16M_CdTe_000364.h5")
dataFile = os.path.join(dataDir,"c103-1-ff-1_*_EIG16M_CdTe_{num2:0>6}.h5")

spotsFile = r"C:\Users\dpqb1\Documents\Data\c103_processing\spots.npz"
spotData = np.load(spotsFile)

params = {}
params['detector'] = 'eiger'
params['peak_func'] = "gaussian"
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,"eiger16M_monolith_mruby_062224_FINAL.yml")
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False

grains = [15,158]
spotInds = sf.findSpots(spotData,grains=grains)
spotInds = [113,205,413,954]
spotInds = [113]

ttPath = os.path.join(topPath,'outputs_test')
# os.mkdir(ttPath)
scanRange = np.concatenate((np.array([364,368,372,376,380]), np.arange(383,406), [407]))
# sf.createTruth(state_file,spotInds[0:2],spotData,scanRange,ttPath,ttPath,dataFile,params)


spotInd = 3
scan0 = 0
print(spotData['etas'][spotInd])
# %% ROI Visualizer
# dome = 3
# num_cols = 8
# sf.roiTrackVisual(spotInd,spotData,dome,num_cols,scanRange,dataFile,ttPath,ttPath,params)
# sf.roiTrackVisual2([spotInd],spotData,dome,scanRange,ttPath,dataFile,params)
# %% Create Truth
root = tk.Tk()
app = sf.truthGUI(root,spotInd,scan0,spotData,scanRange,ttPath,ttPath,dataFile,params)
root.mainloop()