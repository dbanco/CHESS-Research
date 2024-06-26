# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import scipy
import numpy as np
import spotfetch as sf
import chess_detectors as cd
import os
import pickle
import matplotlib.pyplot as plt

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
K = len(spotData['ome_idxs'])
frmArray = spotData['ome_idxs']
dArray = {}

tthLen = 30
etaLen = 100
k = 10

# Load in dexela img stack
paramsD = {}
# params['detector'] = 'eiger'
# params['imSize'] = (5000,5000)
# params['yamlFile'] = '/mnt/scratch/dbanco/c103_processing/eiger16M_monolith_mruby_062224_FINAL.yml'
paramsD['detector'] = 'dexela'
paramsD['imSize'] = (4888,7300) 
paramsD['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
paramsD['roiSize'] = [tthLen,etaLen]
paramsD['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]
paramsD['pool'] = 16

params = {}
params['detector'] = 'eiger'
params['imSize'] = (5000,5000)
params['yamlFile'] = '/mnt/scratch/dbanco/c103_processing/eiger16M_monolith_mruby_062224_FINAL.yml'
# params['detector'] = 'dexela'
# params['imSize'] = (4888,7300) 
# params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
params['roiSize'] = [tthLen,etaLen]
params['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
fnamesE = ['/nfs/chess/id1a3/2024-2/wrangler-orphans/nb-ff-prescans-1_0006_EIG16M_CdTe_000006.h5']

scanNum = 2
fnames = sf.pathToFile(os.path.join(exsituPath,f'{scanNum}','ff'))
# dexImgs = np.zeros((1441,4888,7300))
# roiDexlist = []
# roiEiglist = []
dexList = []
eigList = []
tth = spotData['tths'][k]
eta = spotData['etas'][k]
frm = spotData['ome_dxs'][k]

#%%
dexPolarArray = np.zeros((1441,tthLen,etaLen))
for i in np.arange(4,1445):
    print(i)
    dexPolarArray[i-4,:,:] = cd.loadDexPolarRoi(fnames, tth, eta, i, paramsD)
eigPolar = cd.loadEigerPolarRoi(fnamesE[0],tth,eta,2,params)

outFile = f'dexPolarArraySpot{k}.pkl'
with open(outFile, 'wb') as f:
    pickle.dump(dexPolarArray,f)

# %%
# eigPolar = eigPolar.reshape((1,40,etaLen)
# cc = scipy.ndimage.correlate(dexPolarArray, eigPolar,mode='constant',cval=0)
# print(cc.shape)
# ccMax = np.max(cc,axis=(1,2))
# plt.plot(ccMax)
# omegFrm = np.where(ccMax == np.max(ccMax))[0]
# print(omegFrm)
# print(spotData['ome_idxs'][4])



# %% Extract 2theta, eta region
for i in [k]:
    tth = spotData['tths'][i]
    eta = spotData['etas'][i]
    eigPolar = cd.loadEigerPolarRoi(fnamesE[0],tth,eta,2,params)
    eigPolar = eigPolar.reshape((1,tthLen,100))
    outFile = f'dexPolarArraySpot{i}.pkl'
    with open(outFile,'rb') as f:
        dexPolarArray = pickle.load(f)
      
    cc = np.zeros(1441)
    for i in range(1441):
        dexPolar = dexPolarArray[i,:,:]/np.linalg.norm(dexPolarArray[i,:,:])
        cc[i] = np.sum(np.multiply(eigPolar[0,:,:],dexPolar),axis=None)
    # cc = scipy.ndimage.correlate(dexPolarArray, eigPolar,mode='constant',cval=0)
    
    plt.figure()
    plt.plot(cc)
    omegFrm = np.where(cc == np.max(cc))[0]
    print(omegFrm)
    
# %%
plt.figure()
plt.subplot(2,1,1)
plt.imshow(eigPolar[0,:,:])
plt.title('Eiger polar region')
plt.subplot(2,1,2)
plt.title('Dexela polar region')
plt.imshow(np.squeeze(dexPolarArray[omegFrm,:,:]))

# %% Crop around center to same size
# cropDims = [4888,4888]
# eigImg_crop = eigImg[:,66:]

# %% Cross correlate along dim 0
# max_cc = np.zeros(1441)
# for i in range(1441):
#     max_cc[i] = np.max(scipy.signal.correlate2d(dexImgs[i,:,:], eigImg),axis=None)

# %% c
# for ome_shift in range(1441):
    
#     for k in range(K):
#         frmArray[k] = sf.wrapFrame(frmArray[k] - 235 ,frm0=4,numFrms=1441)
#     dArray['ome_idxs'] = frmArray
#     dArray['etas'] = spotData['etas']
#     dArray['tths'] = spotData['tths']
#     # %% 3. Detector and tracking parameters
#     params = {}
#     params['detector'] = 'eiger'
#     params['imSize'] = (5000,5000)
#     params['yamlFile'] = '/mnt/scratch/dbanco/c103_processing/eiger16M_monolith_mruby_062224_FINAL.yml'
#     # params['detector'] = 'dexela'
#     # params['imSize'] = (4888,7300) 
#     # params['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
#     params['roiSize'] = 40
#     params['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]
#     params['pool'] = 16
    
#     # %% 4. Inspect spot tracks on ex-situ and initial scan data
#     scanNum = 2
#     fNum = 96
#     frame = 4
#     detectFrame = 2
#     # fnames = sf.pathToFile(os.path.join(exsituPath,f'{scanNum}','ff'))
#     fnames = ['/nfs/chess/id1a3/2024-2/wrangler-orphans/nb-ff-prescans-1_0006_EIG16M_CdTe_000006.h5']
#     sf.plotSpotWedges(dArray,fnames,frame,params,detectFrame)
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