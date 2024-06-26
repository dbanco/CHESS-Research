# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os
import matplotlib.pyplot as plt
import chess_detectors as cd
import scipy
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
paramsE = {}
paramsE['detector'] = 'eiger'
paramsE['imSize'] = (5000,5000)
paramsE['yamlFile'] = '/mnt/scratch/dbanco/c103_processing/eiger16M_monolith_mruby_062224_FINAL.yml'
paramsE['roiSize'] = [40,40]
paramsE['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]
paramsE['pool'] = 16

paramsD = {}
paramsD['detector'] = 'dexela'
paramsD['imSize'] = (4888,7300) 
paramsD['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
paramsD['roiSize'] = [40,40]
paramsD['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]
paramsD['pool'] = 16

# %% 4. Inspect spot tracks on ex-situ and initial scan data
scanNum = 2
fNum = 96
dexFrame = 4
eigFrames = np.arange(2,1442) # Checked 224 to240
fnamesD = sf.pathToFile(os.path.join(exsituPath,f'{scanNum}','ff'))
fnamesE = ['/nfs/chess/id1a3/2024-2/wrangler-orphans/nb-ff-prescans-1_0006_EIG16M_CdTe_000006.h5']
# sf.plotSpotWedges(spotData,fnamesD,dexFrame,paramsD,dexFrame)
# for eigFrame in eigFrames:
#     sf.plotSpotWedges(spotData,fnamesE,dexFrame,paramsE,eigFrame)
#     plt.title(f'Eiger frame {eigFrame}')
    
dexImg = cd.load_dex(fnamesD,paramsD,dexFrame)
# eigImgArray = np.zeros((1440,5000,5000))
# for eigFrame in eigFrames:
#     print(eigFrame)
#     eigImgArray[eigFrame-2,:,:] = cd.load_eiger(fnamesE,paramsE,eigFrame)


# %% Pick out the side of raing region
# np.load('eigImgArray.npy')
dexSubImg = dexImg[1600:3200,1600:2400]
# eigSubArray = eigImgArray[:,1900:3100,1000:1600]

# plt.figure()
# plt.subplot(1,2,1)
# plt.title('Eiger {eigFrm}')
# eigSubImg = np.squeeze(eigSubArray[232,:,:])
# plt.imshow(eigSubImg)
# plt.clim(np.median(eigSubImg),np.median(eigSubImg)+100)
# plt.colorbar()

# plt.subplot(1,2,2)
# plt.title('Dexela {dexFrame}')
# plt.imshow(dexSubImg)
# plt.clim(np.median(dexSubImg),np.median(dexSubImg)+100)
# plt.colorbar()

# %% Concatente a bunch of eiger imaegs

dexSubImg = dexImg[1600:3200,1600:2400]

frm1 = 751
frm2 = 752
eigImg = cd.load_eiger(fnamesE,paramsE,frm1)
# eigConcat = eigImg[1900:3100,1000:1600]
for frmi in np.arange(frm1,frm2):
    eigImg = cd.load_eiger(fnamesE,paramsE,frmi)
    # eigConcat = np.hstack((eigConcat,eigImg[1900:3100,1000:1600]))
    
    plt.figure()
    plt.title(f'Eiger {frmi}')
    plt.imshow(eigImg[1900:3100,1000:1600])
    plt.clim(np.median(eigImg[1900:3100,1000:1600]),np.median(eigImg[1900:3100,1000:1600])+100)
    plt.colorbar()

plt.figure()
plt.title(f'Dexela {dexFrame}')
plt.imshow(dexSubImg)
plt.clim(np.median(dexSubImg),np.median(dexSubImg)+100)
plt.colorbar()

# Checked 228-236 (56.5-58.5)
# Checekd 588-595 (146.5-148.5)
# (maybe 595)
# Checked 748-756 (186.5-188.5)
# (maybe 751)


# %%
# eigSubArray = eigImgArray[:,1900:3100,1000:1600]
# cc = scipy.ndimage.correlate(eigSubArray, dexSubImg.reshape((1,1600,800)),mode='constant',cval=0)
# ccMax = np.max(cc,axis=(1,2))
# plt.figure()
# plt.plot(ccMax)
# omegFrm = np.where(ccMax == np.max(ccMax))[0]
# print(omegFrm)
# print(spotData['ome_idxs'][k])


# %%

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