# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os
# import multiprocesing

#   1. Set up paths: exsitu, data, outputs
topPath = "/mnt/scratch/dbanco/c103_processing/Sample-1/layer-1"
exsituPath = "/nfs/chess/raw/2024-2/id3a/miller-3528-c/c103-1-ff-1/"
dataPath = "/mnt/scratch/nygren-series-cycle-2024-2/"

# %% 2. Load in or collect spots data
# spotsDir = "TBD"
spotsPath = os.path.join(topPath,'spots')
# Get data from spots data
# sf.collectSpotsData(spotsPath, spotsDir) 
spotData = np.load(os.path.join(spotsPath,'spots.npz'))
# K = len(spotData['ome_idxs'])
# frmArray = spotData['ome_idxs']
# dArray = {}
# for k in range(K):
#     frmArray[k] = sf.wrapFrame(frmArray[k] + 557 ,frm0=4,numFrms=1441)
# dArray['ome_idxs'] = frmArray
# dArray['etas'] = spotData['etas']
# dArray['tths'] = spotData['tths']
# %% 3. Detector and tracking parameters
params = {}
params['detector'] = 'eiger'
params['imSize'] = (5000,5000)
params['yamlFile'] = '/mnt/scratch/dbanco/c103_processing/eiger16M_monolith_mruby_062224_FINAL.yml'
# params['detector'] = 'dexela'
# params['imSize'] = (4888,7300) 
# params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
params['roiSize'] = [40,40]
params['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16

# %% 4. Inspect spot tracks on ex-situ and initial scan data
scanNum = 2
fNum = 96
frame = 4
detectFrame = 1132
# fnames = sf.pathToFile(os.path.join(exsituPath,f'{scanNum}','ff'))
fnames = ['/nfs/chess/id1a3/2024-2/wrangler-orphans/nb-ff-prescans-1_0006_EIG16M_CdTe_000006.h5']
sf.plotSpotWedges(spotData,fnames,frame,params,detectFrame)
# roi_list = sf.loadSpotsAtFrame(dArray,fnames,frame,params,detectFrame)
# sf.plotROIs(roi_list)

# scanNum = 1
# fNum = 95
# frame = 10
# fname1 = os.path.join(dataPath,f'{scanNum}','ff','ff1_000' + str(fNum) + '.h5')
# fname2 = os.path.join(dataPath,f'{scanNum}','ff','ff2_000' + str(fNum) + '.h5')
# sf.plotSpotWedges(spotData,fname1,fname2,frame,params)

# # %% 5. Determine which spots to process
# # spotInds = np.arange(120,1000)
# grains = [5,11]
# tth_ring = 3.76*np.pi/180
# dtth_ring= 4*np.pi/180
# spotInds = sf.findSpots(spotData,grains=grains,tth=tth_ring,dtth=dtth_ring)

# # %% 6. Begin Processing
# scan1 = 1
# initTracksPath = os.path.join(topPath,'outputs')
# sf.initExsituTracks(initTracksPath,exsituPath, spotData, spotInds, params, scan1-1)
# sf.spotTracker(dataPath,topPath,exsituPath,spotData,spotInds,params,scan1)