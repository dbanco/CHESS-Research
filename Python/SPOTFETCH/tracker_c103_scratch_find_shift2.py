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
from scipy import stats

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

tthLen = 40
etaLen = 40
ccList = []
for k in range(20):
    print(k)
    tth = spotData['tths'][k]
    eta = spotData['etas'][k]
    frm = spotData['ome_idxs'][k]
    
    
    # Load in dexela img stack
    paramsD = {}
    paramsD['detector'] = 'dexela'
    paramsD['imSize'] = (4888,7300) 
    paramsD['yamlFile'] = '/nfs/chess/user/dbanco/c103_processing/dexelas_calibrated_ruby_0504_v01.yml'
    paramsD['roiSize'] = [tthLen,etaLen]
    paramsD['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]
    paramsD['pool'] = 16
    scanNum = 2
    # fnames = sf.pathToFile(os.path.join(exsituPath,f'{scanNum}','ff'))
    # dex_img = cd.load_dex(fnames,paramsD,frm)
    # plt.imshow(dex_img)
    # plt.clim(300,400)
    
    # %% 
    params = {}
    params['detector'] = 'eiger'
    params['imSize'] = (5000,5000)
    params['yamlFile'] = '/mnt/scratch/dbanco/c103_processing/eiger16M_monolith_mruby_062224_FINAL.yml'
    params['roiSize'] = [tthLen,etaLen]
    params['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]
    params['pool'] = 16
    fnamesE = ['/nfs/chess/id1a3/2024-2/wrangler-orphans/nb-ff-prescans-1_0006_EIG16M_CdTe_000006.h5']
    
    # detectDistE, mmPerPixelE, ff_transE = cd.loadYamlData(params, tth, eta)
    # detectDistD, mmPerPixelD, ff_transD = cd.loadYamlData(paramsD, tth, eta)
    
    # # cd.polarDomain(detectDist,mmPerPixel,tth,eta,roi_size)
    # rE = np.round(detectDistE*np.tan(tth)/mmPerPixelE)
    # r2E = rE + tthLen//2
    # detaE = 1/r2E
    # etaWidE = detaE*etaLen
    
    # rD = np.round(detectDistD*np.tan(tth)/mmPerPixelD)
    # r2D = rD + tthLen//2
    # detaD = 1/r2D
    # etaWidD = detaD*etaLen
    
    # eta1 = eta - roi_size[1]/2*deta
    # eta2 = eta + roi_size[1]/2*deta
    # eta_domain = np.linspace(eta1,eta2,roi_size[1])
    
    # %% 4. Inspect spot tracks on ex-situ and initial scan data
    scanNum = 2
    fNum = 96
    eigFrame = 1117 # spot 37 gives 1117
    detectFrame = 2
    # fnames = sf.pathToFile(os.path.join(exsituPath,f'{scanNum}','ff'))
    # fnames = ['/nfs/chess/id1a3/2024-2/wrangler-orphans/nb-ff-prescans-1_0006_EIG16M_CdTe_000006.h5']
    # sf.plotSpotWedges(dArray,fnames,frame,params,detectFrame)
    
    scanNum = 2
    fnames = sf.pathToFile(os.path.join(exsituPath,f'{scanNum}','ff'))
    
    dexPolar = cd.loadDexPolarRoi(fnames,tth,eta,frm,paramsD)
    eigPolar = cd.loadEigerPolarRoi(fnamesE[0], tth, eta, eigFrame, params)
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.imshow(dexPolar)
    plt.colorbar()
    plt.title('Dexela polar region')
    plt.clim(np.median(dexPolar),np.median(dexPolar)+100)
    plt.subplot(2,1,2)
    plt.imshow(eigPolar)
    plt.colorbar()
    plt.title('Eiger polar region')
    plt.clim(1,20)
    
    # %% SUm in tth
    # eigOmgEta = np.sum(eigPolarArray,1)
    # plt.figure()
    # plt.imshow(np.transpose(eigOmgEta))
    
    #%%
    eigPolarArray = cd.loadEigerPolarRoiArray(fnamesE[0], tth, eta, [2,1442], params)
    # np.load('eigPolarArray.npy')
    # outFile = f'eigPolarArraySpot{k}.pkl'
    # with open(outFile, 'wb') as f:
    #     pickle.dump(eigPolarArray,f)
    for i in range(eigPolarArray.shape[0]):
        eigPolarArray[i,:,:] = eigPolarArray[i,:,:]/np.linalg.norm(eigPolarArray[i,:,:])
    dexPolarRe = dexPolar.reshape((1,tthLen,etaLen))
    cc = scipy.ndimage.correlate(eigPolarArray, dexPolarRe,mode='constant',cval=0)
    ccList.append(cc)
    
# %%
ccSum = np.zeros(cc.shape)
omegArray = np.zeros(20)
for i,cc_i in enumerate(ccList):
    ccSum += cc_i
    ccMax_i = np.max(cc_i,axis=(1,2))
    ccMax2 = np.max(ccMax_i)
    omegFrm = np.where(ccMax_i == ccMax2)[0]
    omegArray[i] = omegFrm
    
ccMax = np.max(ccSum,axis=(1,2))
omegFrmSum = np.where(ccMax == np.max(ccMax))[0]

plt.figure()
plt.plot(omegArray)

plt.figure()
plt.plot(ccMax)

print(omegFrmSum)
print(spotData['ome_idxs'][k])

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
# for i in [k]:
#     tth = spotData['tths'][i]
#     eta = spotData['etas'][i]
#     frm = spotData['ome_idxs'][i]
#     dexPolar = cd.loadDexPolarRoi(fnames,tth,eta,frm,paramsD)
#     outFile = f'eigPolarArraySpot{i}.pkl'
#     with open(outFile,'rb') as f:
#         eigPolarArray = pickle.load(f)
      
#     cc = np.zeros(1440)
#     for i in range(1440):
#         eigPolar = eigPolarArray[i,:,:]/np.linalg.norm(eigPolarArray[i,:,:])
#         cc[i] = np.sum(np.multiply(eigPolar,dexPolar),axis=None)
#     # cc = scipy.ndimage.correlate(dexPolarArray, eigPolar,mode='constant',cval=0)
    
#     plt.figure()
#     plt.plot(cc)
#     omegFrm = np.where(cc == np.max(cc))[0]
#     print(omegFrm)
    
# # %%
# plt.figure()
# plt.subplot(2,1,1)
# plt.imshow(np.squeeze(eigPolarArray[783,:,:]))
# plt.title('Eiger polar region')
# plt.subplot(2,1,2)
# plt.title('Dexela polar region')
# plt.imshow(dexPolar)


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