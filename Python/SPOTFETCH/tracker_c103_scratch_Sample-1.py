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
topPath = "/mnt/scratch/dbanco/c103_processing/Sample-1"
# exsituPath = "/nfs/chess/raw/2024-2/id3a/miller-3528-c/c103-1-ff-1/"
dataPath = "/mnt/scratch/nygren-series-cycle-2024-2/"
# dataPath = "/nfs/chess/id1a3/2024-2/nygren-4125-a/nygren-series-cycle-2024-2-chessdaq/nygren-series-cycle-2024-2"
num1 = 2
num2 = 353
# dataFile = os.path.join(dataPath,f'c103-2-ff-1_{num1:0>{4}}_EIG16M_CdTe_{num2:0>{6}}.h5')
# exsituPath = os.path.join(dataPath,f'c103-3-ff-1_{num1:0>{4}}_EIG16M_CdTe_{num2:0>{6}}.h5')
exsituPath = os.path.join(dataPath,f'c103-1-ungripped-1_{num1:0>4}_EIG16M_CdTe_{num2:0>6}.h5')
dataFile = os.path.join(dataPath,"c103-1-ff-1_*_EIG16M_CdTe_{num2:0>6}.h5")
# dataFile = os.path.join(dataPath,"c103-1-ff-1_*_EIG16M_CdTe_{num2:0>6}.h5")


# %% 2. Load in or collect spots data
spotsDir = '/nfs/chess/user/wkk32/c103-1-layer1/'
 
spotsPath = os.path.join(topPath,'spots')
# Get data from spots data
sf.collectSpotsData(spotsPath, spotsDir) 
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
params['pool'] = 100

# %% 4. Inspect spot tracks on ex-situ and initial scan data
grains = [276,288,342]
spotInds = sf.findSpots(spotData,grains=grains)

dome = 3
scanRange = np.concatenate((np.array([353, 364,368,372,376,380]), np.arange(383,406), [407]))
trackPath = os.path.join(topPath,'outputs')

sf.roiTrackVisual(np.arange(383,406),spotData,dome,scanRange,trackPath,dataFile,params)

frame = 17
# sf.plotSpotWedges(spotData,exsituPath,frame,params,grains=grains)

# roi_list = sf.loadSpotsAtFrame(dArray,fnames,frame,params,detectFrame)
# sf.plotROIs(roi_list)s

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

# tthdeg = 4
# tth_ring = tthdeg*np.pi/180
# dtth_ring= 0.5*np.pi/180
# eta1 = np.pi/2
# eta2 = 0
# grains1 = [5]
# spotInds1 = sf.findSpots(spotData,eta=np.pi/2,deta=0.25,tth=tth_ring,dtth=tth_ring)
# spotInds2 = sf.findSpots(spotData,eta=np.pi,deta=10,tth=tth_ring,dtth=dtth_ring)

# # %% 6. Begin Processing
initTracksPath = os.path.join(topPath,'outputs')
# sf.initExsituTracks(initTracksPath,exsituPath, spotData, spotInds, params, 352)

# Sequence: 353, 364, 368, 372, 376, 380
num1 = 4
num2 = 383
advance = True
sf.spotTracker(dataFile,topPath,spotData,spotInds,params,num1,num2,advance)