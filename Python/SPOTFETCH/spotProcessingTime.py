# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 09:29:21 2024

@author: dpqb1
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import spotfetch as sf
import processingSigmaXY as pXY
from hexrd.fitting import fitpeak
from hexrd.fitting import peakfunctions as pkfuncs
import time

# %% First load indexed spot data
folder_path = "spots_11032023"  
grain_data = np.load(folder_path+'.npz')    
Xs = grain_data['Xs']
Ys = grain_data['Ys']
id_nums = grain_data['id_nums']
tths = grain_data['tths']
etas = grain_data['etas']    
omes = grain_data['omes']    
ome_idxs = grain_data['ome_idxs']  
grain_nums = grain_data['grain_nums']

# %%  Load data for 4th frame for all timesteps
detectDist = 608
mmPerPixel = 0.0748
center = [4888/2,7300/2]


fDir = "D:\\CHESS_raw_data\\ti-2-tension\\"
fPath = fDir+'12\\ff'
fname1 = fDir + '12\\ff\\ff1_000098.h5'
fname2 = fDir + '12\\ff\\ff2_000098.h5'
fileNum = 106
numDirs = 73

frame = 4
spotInds = np.where(ome_idxs==frame)[0]
num_rois = len(spotInds)
roi_size = 40
fit_params_hexrd = np.zeros([num_rois,6,numDirs])
fit_params_lmfit = np.zeros([num_rois,6,numDirs])

# %%  Run this block to continue processing without reseting
for dirNum in range(61,numDirs):
    dNum = str(dirNum + 1)
    fNum = str(fileNum)
    fname1 = fDir + dNum + '\\ff\\ff1_000' + fNum + '.h5'
    fname2 = fDir + dNum + '\\ff\\ff2_000' + fNum + '.h5'
    print(dNum)
    
    # Load full HEXD image  
    try:
        b = sf.load_dex(fname1,fname2,frame-1,frame+1)
    except:
        continue

    # Load ROIs
    rois = np.zeros((roi_size,roi_size,num_rois))
    j = 0
    for i in spotInds:
        tth = tths[i]
        eta = etas[i]
        rois[:,:,j] = sf.get_roi(b,detectDist,mmPerPixel,center,tth,eta,roi_size)
        # sf.add_spot_rect(detectDist,mmPerPixel,center,tth,eta,roi_size)
        j += 1
    
    # Process Spots
    for i in range(num_rois):      
        # ROI image  (Could preprocess full image instead of patches)
        img = rois[:,:,i]
        img_median = np.median(img)
        img = img - 1.1*img_median
        img[img < 0] = 0

        # HEXRD Peak fitting
        # Uses scipy.optimize.leastsq
        # Also estimates the background parameters. This cannot be removed
        # (bg0+bg1x*x+bg1y*y)
        y_vals, x_vals = np.indices(img.shape)
        p0 = fitpeak.estimate_pk_parms_2d(x_vals, y_vals, img,'gaussian')
        p0[5:8]=0
        p = fitpeak.fit_pk_parms_2d(p0, x_vals, y_vals, img, 'gaussian')
        if p[3]==0: p[3] += 1e-4
        if p[4]==0: p[4] += 1e-4
        fit_params_hexrd[i,0:5,dirNum] = p[0:5]
        recon = pkfuncs.gaussian2d(p, x_vals, y_vals)
        fit_params_hexrd[i,5,dirNum] = np.linalg.norm(recon-img)/np.linalg.norm(img)
        
        # LMfit Peak fitting
        x, y = pXY.findBlob(roi_size//2, roi_size//2, img)
        result = pXY.fitGaussian(img, x, y)
        fit_params_lmfit[i,0:5,dirNum] = pXY.getGaussFitParams(result)
        recon = pXY.gaussianRecon(result, roi_size)
        fit_params_lmfit[i,5,dirNum] = np.linalg.norm(recon-img)/np.linalg.norm(img)
           
    fileNum += 1

# %% Plot FWHM time series
plt.figure(6)
for i in range(num_rois):
    plt.subplot(5,5,i+1)
    plt.plot(fit_params_hexrd[i,3,:],marker='x',color='r',markersize=10)
    plt.plot(fit_params_lmfit[i,3,:],marker='x',color='c',markersize=10)
    plt.xlabel('Time Step')
    plt.ylabel('FWHM x')
    
# %% View particular ROI
plt.figure(7)
k = 1
frames = [51,52,53,54,55,56,57,58]
for timeStep in frames:
    roiNum = 19
    
    dNum = str(timeStep + 1)
    fNum = str(106 + timeStep)
    fname1 = fDir + dNum + '\\ff\\ff1_000' + fNum + '.h5'
    fname2 = fDir + dNum + '\\ff\\ff2_000' + fNum + '.h5'
    
    # Load full HEXD image  
    b = sf.load_dex(fname1,fname2,frame-1,frame+1)

    # Load ROIs
    rois = np.zeros((roi_size,roi_size,num_rois))
    j = 0
    for i in spotInds:
        tth = tths[i]
        eta = etas[i]
        rois[:,:,j] = sf.get_roi(b,detectDist,mmPerPixel,center,tth,eta,roi_size)
        # sf.add_spot_rect(detectDist,mmPerPixel,center,tth,eta,roi_size)
        j += 1

    img = rois[:,:,roiNum]
    # img_median = np.median(img)
    # img = img - 1.1*img_median
    # img[img < 0] = 0

    plt.subplot(1,len(frames),k)
    plt.imshow(img)
    k += 1