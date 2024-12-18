# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os
import glob
import matplotlib.pyplot as plt

# LOCAL
topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
exsituPath = os.path.join(topPath,r"state_0\simulation\outputs\c103_polycrystal_sample_2_state_0_layer_1_output_data.npz")
dataFile = os.path.join(topPath,r"state_*\simulation\outputs\c103_polycrystal_sample_2_state_*_layer_1_output_data.npz")

# %% 2. Load in or collect spots data
# Spots for each state file
state = 4
spotsDir = os.path.join(topPath,f'state_{state}','simulation','outputs')
spotsOut = os.path.join(topPath,f'state_{state}')

# sf.collectSpotsData(spotsOut, spotsDir) 
spotData = np.load(os.path.join(spotsOut,'spots.npz'))

# %% Detector and tracking parameters
params = {}
params['detector'] = 'eiger_sim'
params['peak_func'] = "gaussian"
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,'c103_eiger_calibration.yml')
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False

fullscanRange = np.arange(5)
trackPath = os.path.join(topPath,'outputs_test')
fname = os.path.join(topPath,f"state_{state}\simulation\outputs\c103_polycrystal_sample_2_state_{state}_layer_1_output_data.npz")

# sf.plotSpotWedges(spotData,fname,spotData['ome_idxs'][0],params)

# spotInd = 0
# dome = 2
# sf.roiTrackVisual(spotInds,spotData,dome,fullscanRange,trackPath,dataFile,params)
# sf.roiTrackVisual(spotInd,spotData,dome,fullscanRange,dataFile,trackPath,trackPath,params)

Nrows = 5
Ncols = 5

# fig, axes = plt.subplots(Nrows,Ncols, figsize=(10, 10))
# k = 300
# for ax in axes.ravel():
#     tthRoi = spotData['tths'][k]
#     etaRoi = spotData['etas'][k]
#     sf.showROI(ax,fname,0,
#                 spotData['ome_idxs'][k],
#                 tthRoi,
#                 etaRoi,params)
#     sf.showInitial(ax,etaRoi,tthRoi,etaRoi,tthRoi,params)
#     k = k + 1
    

fig, axes = plt.subplots(Nrows,Ncols, figsize=(10, 10))
k = 300
for ax in axes.ravel():
    x_meas = spotData['Xm'][k]
    y_meas = spotData['Ym'][k]
    eta,tth = sf.xyToEtaTthRecenter(x_meas, y_meas, params)
    sf.showROI(ax,fname,0,spotData['ome_idxs'][k],tth,eta,params)
    sf.showInitial(ax,eta,tth,eta,tth,params)
    k = k + 1
    