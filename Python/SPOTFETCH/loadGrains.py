# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:44:45 2024

@author: dpqb1
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import spotfetch as sf
import processingSigmaXY as pXY
from hexrd.fitting import fitpeak
from hexrd.fitting import peakfunctions as pkfuncs
import time
# from hexrd import imageseries
# from hexrd.imageseries.process import ProcessedImageSeries as ProcessedIS
# from hexrd.imageseries.omega import OmegaWedges


# %% Define data file and indexed grain file
fDir = "D:\\CHESS_raw_data\\ti-2-exsitu\\"
fPath = fDir+'12\\ff'
fname1 = fDir + '12\\ff\\ff1_000098.h5'
fname2 = fDir + '12\\ff\\ff2_000098.h5'

# fDir = "D:\\CHESS_raw_data\\ti-2-tension\\"
# fPath = fDir+'65\\ff'
# fname1 = fDir + '65\\ff\\ff1_000169.h5'
# fname2 = fDir + '65\\ff\\ff2_000169.h5'


folder_path = "spots_11032023"  
# folder_path = "spots_files"
# folder_path = "spots_11022023"
# folder_path = "spots_11132023"  

sf.collectSpotsData(folder_path)

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

yamlFile = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\dex-refined-1.yml"
detectDist = 608
mmPerPixel = 0.0748
center = [4888/2,7300/2]


# %% Load stack of ROIs for first omega
# frame = 4
# b = sf.load_dex(fname1,fname2,frame-1,frame+1)
# # plt.figure()
# # plt.imshow(b)
# # plt.plot(center[1],center[0],marker="x",color="y")
# # plt.clim(290,550)

# num_rois = len(np.where(ome_idxs==frame)[0])
# roi_size = 40
# rois = np.zeros((roi_size,roi_size,num_rois))
# j = 0
# for i in np.where(ome_idxs==frame)[0]:
#     tth = tths[i]
#     eta = etas[i]
#     rois[:,:,j] = sf.get_roi(b,detectDist,mmPerPixel,center,tth,eta,roi_size)
#     # sf.add_spot_rect(detectDist,mmPerPixel,center,tth,eta,roi_size)
#     j += 1

    
# %% Fitting with HEXRD Gaussian code and Connor's Code
# fit_params_hexrd = np.zeros([num_rois,5])
# fit_params_connor = np.zeros([num_rois,5])
# for i in range(num_rois):
    
#     # ROI image  (Could preprocess full image instead of patches)
#     img = rois[:,:,i]
#     img_median = np.median(img)
#     img = img - 1.1*img_median
#     img[img < 0] = 0
#     plt.figure(1)
#     plt.subplot(5,5,i+1)
#     plt.imshow(img)
        
#     # HEXRD Peak fitting
#     # Uses scipy.optimize.leastsq
#     # Also estimates the background parameters. This cannot be removed
#     # (bg0+bg1x*x+bg1y*y)
#     avg_time1= 0
#     tic = time.perf_counter()
#     y_vals, x_vals = np.indices(img.shape)
#     p0 = fitpeak.estimate_pk_parms_2d(x_vals, y_vals, img,'gaussian')
#     p0[5:8]=0
#     p = fitpeak.fit_pk_parms_2d(p0, x_vals, y_vals, img, 'gaussian')
#     if p[3]==0: p[3] += 1e-4
#     if p[4]==0: p[4] += 1e-4
#     fit_params_hexrd[i] = p[0:5]
    
#     recon = pkfuncs.gaussian2d(p, x_vals, y_vals)
#     toc = time.perf_counter()
#     avg_time1 += (toc-tic)/num_rois
#     plt.figure(2)
#     plt.subplot(5,5,i+1)
#     plt.imshow(recon)
    
#     # Connor Peak fitting
#     avg_time2 = 0
#     tic = time.perf_counter()
#     x, y = pXY.findBlob(roi_size//2, roi_size//2, img)
#     result = pXY.fitGaussian(img, x, y)
#     fit_params_connor[i] = pXY.getGaussFitParams(result)
#     recon = pXY.gaussianRecon(result, roi_size)
#     toc = time.perf_counter()
#     avg_time2 += (toc-tic)/num_rois
#     plt.figure(3)
#     plt.subplot(5,5,i+1)
#     plt.imshow(recon)
    
#     plt.figure(1)
#     plt.subplot(5,5,i+1)
#     plt.plot(fit_params_hexrd[i,1],fit_params_hexrd[i,2],marker='x',color='r',markersize=10)
#     plt.plot(fit_params_connor[i,1],fit_params_connor[i,2],marker='x',color='c',markersize=10)

# print(f"HEXRD: {avg_time1:0.4f} seconds")
# print(f"Connor: {avg_time2:0.4f} seconds")\

# %% Plot FWHM pairs spot by spot
# plt.figure(4)
# plt.subplot(2,1,1)
# plt.plot(fit_params_hexrd[:,3],marker='x',color='r')
# plt.plot(fit_params_connor[:,3],marker='o',color='c',fillstyle='none')
# plt.xlabel('Spot Index')
# plt.ylabel('FWHM_x')
# plt.legend('HEXRD','LMfit')

# plt.subplot(2,1,2)
# plt.plot(fit_params_hexrd[:,4],marker='x',color='r')
# plt.plot(fit_params_connor[:,4],marker='o',color='c',fillstyle='none')
# plt.xlabel('Spot Index')
# plt.ylabel('FWHM_y')
# plt.legend('HEXRD','LMfit')

# Plot FWHM for single spot over time


# %% Plot rectangles for grain 0
'''
b = sf.load_dex(fname1,fname2,1345)

grain0 = np.where(grain_nums==0)[0]
# omeSel = np.where((ome_idxs==1344))[0]
grain0omeSel = np.where((ome_idxs==1346) & (grain_nums==0))[0]

plt.figure()
plt.imshow(b)
plt.clim(260,550)
plt.plot(center[1],center[0],marker="x",color="y")

for j in grain0omeSel:
    tth = tths[j]
    eta = etas[j]
    print(ome_idxs[j])
    sf.add_spot_rect(detectDist,mmPerPixel,center,tth,eta)
'''