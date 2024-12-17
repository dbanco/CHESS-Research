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

topPath = r"C:\Users\dpqb1\Documents\Data\c103_processing"
dataDir = r"E:\Data\c103"

exsituPath = os.path.join(dataDir,"c103-1-ff-1_0004_EIG16M_CdTe_000364.h5")
dataFile = os.path.join(dataDir,"c103-1-ff-1_*_EIG16M_CdTe_{num2:0>6}.h5")

spotsFile = r"C:\Users\dpqb1\Documents\Data\c103_processing\spots.npz"
spotData = np.load(spotsFile)

# Detector and tracking parameters
params = {}
params['detector'] = 'eiger'
params['peak_func'] = "gaussian"
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,"eiger16M_monolith_mruby_062224_FINAL.yml")
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False

fullscanRange = np.concatenate((np.array([364,368,372,376,380]), np.arange(383,406),[407]))
trackPath = os.path.join(topPath,'outputs_test')

template = dataFile.format(num2=364)
fname = glob.glob(template)[0]

# sf.plotSpotWedges(spotData,fname,spotData['ome_idxs'][0],params)


# spotInd = 0
# dome = 2
# sf.roiTrackVisual(spotInds,spotData,dome,fullscanRange,trackPath,dataFile,params)
# sf.roiTrackVisual(spotInd,spotData,dome,fullscanRange,dataFile,trackPath,trackPath,params)

spotInd = 4
grain_id = spotData['id_nums'][spotInd]
spot35 = np.where(spotData['id_nums'] == grain_id)[0]

spotData['omes'][4]
spotData['ome_idxs'][4]


fig, axes = plt.subplots(5,5, figsize=(10, 10))
i = 0
for ax in axes.ravel():
    k = spot35[i]
    x_meas = spotData['Xm'][k]
    y_meas = spotData['Ym'][k]
    eta,tth = sf.xyToEtaTthRecenter(x_meas, y_meas, params)
    sf.showROI(ax,fname,364,spotData['ome_idxs'][k],tth,eta,params)
    i = i + 1