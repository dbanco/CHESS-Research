# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os

topPath = r"C:\Users\dpqb1\Documents\Data\c103_processing"
dataDir = r"D:\Data\c103"

exsituPath = os.path.join(dataDir,"c103-1-ff-1_0004_EIG16M_CdTe_000364.h5")
dataFile = os.path.join(dataDir,"c103-1-ff-1_*_EIG16M_CdTe_{num2:0>6}.h5")

spotsFile = r"C:\Users\dpqb1\Documents\Data\c103_processing\spots.npz"
spotData = np.load(spotsFile)

# Detector and tracking parameters
params = {}
params['detector'] = 'eiger'
params['peak_func'] = 'Gaussian'
params['imSize'] = (5000,5000)
params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/eiger16M_monolith_mruby_062224_FINAL.yml'
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]0
params['pool'] = 16

# Inspect spot tracks on ex-situ and initial scan data
fullscanRange = np.concatenate((np.array([364,368,372,376,380]), np.arange(383,408)))
trackPath = os.path.join(topPath,'outputs')

spotInds = [953]
dome = 3
sf.roiTrackVisual(spotInds,spotData,dome,fullscanRange,trackPath,dataFile,params)