# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os

#   1. Set up paths: exsitu, data, outputs
topPath = "/nfs/chess/user/dbanco/c103_processing/Sample-1"
dataPath = "/nfs/chess/id1a3/2024-2/nygren-4125-a/nygren-series-cycle-2024-2-chessdaq"
num1 = 2
num2 = 353
# exsituPath = os.path.join(dataPath,f'c103-1-ungripped-1_{num1:0>4}_EIG16M_CdTe_000353.h5')
exsituPath = os.path.join(dataPath,"c103-1-ff-1_0004_EIG16M_CdTe_000364.h5")
dataFile = os.path.join(dataPath,"c103-1-ff-1_*_EIG16M_CdTe_{num2:0>6}.h5")


# %% 2. Load in or collect spots data
spotsDir = '/nfs/chess/user/wkk32/c103-1-layer1/'
spotsPath = os.path.join(topPath,'spots')

# The new spots Wiley shared are at
spotsDir = '/nfs/chess/user/dbanco/C103_1_unloaded_gripped_layer2of4/layer_1/3A_grains/'
spotsPath = os.path.join(spotsDir,'eiger')

#sf.collectSpotsData(spotsDir, spotsDir) 
spotData = np.load(os.path.join(spotsDir,'spots.npz'))

# %% 3. Detector and tracking parameters
params = {}
params['detector'] = 'eiger'
params['peak_func'] = 'Gaussian'
params['imSize'] = (5000,5000)
params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/eiger16M_monolith_mruby_062224_FINAL.yml'
# params['detector'] = 'dexela'
# params['imSize'] = (4888,7300) 
# params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]0
params['pool'] = 16

# %% 4. Inspect spot tracks on ex-situ and initial scan data
grains = [44,158]
tth110 = 0.06562438
tth200 = 0.09267698
tth211 = 0.1136209
tth321 = 0.1736603
tths = [tth110, tth211]
hklNames = ['110', '211']
spotInds = sf.findSpots(spotData,grains=grains, tth=tths, dtth=0.012)
# spotInds = [113,205,413,801]

dome = 3
scanRange = np.concatenate((np.array([364,368,372,376,380]), np.arange(383,406), [407]))
trackPath = os.path.join(topPath,'outputs_9_26')

# frame = 17
# sf.plotSpotWedges(spotData,exsituPath,frame,params,grains=grains)

# roi_list = sf.loadSpotsAtFrame(dArray,fnames,frame,params,detectFrame)
# sf.plotROIs(roi_list)s

# # %% 6. Begin Processing
initTracksPath = os.path.join(topPath,'outputs_9_26')
sf.initExsituTracks(initTracksPath,exsituPath,spotData, spotInds, params, 364)

advance = False
# scanRange = np.concatenate((np.array([368,372,376,380]), np.arange(383,406), [407]))
scanRange = np.concatenate(np.arange(397,406), [407]))
for num2 in scanRange:
    sf.spotTracker(dataFile,topPath,spotData,spotInds,params,num1,num2,advance)

# spotInds = [113,205,413,801]
# sf.roiTrackVisual(spotInds,spotData,dome,scanRange,trackPath,dataFile,params)