# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:20:06 2024

@author: dpqb1
"""

import numpy as np
import matplotlib.pyplot as plt
from hexrd.fitting import fitpeak
# import os
import spotfetch as sf
# import processingSigmaXY as pXY
# from hexrd.fitting import fitpeak
# from hexrd.fitting import peakfunctions as pkfuncs
# import time
# import h5py

# %% Test of loading data that  Marianne shared
# aroiDir = 'C:\\Users\\dpqb1\\Documents\\Data\\AROI_data\\Image Series 53'

# f = h5py.File(os.path.join(aroiDir,'Image_series_53_AROI_image_0.mat'),'r')
# img = f['AROI_data']

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

# %% Dataset Parameters
fDir = "D:\\CHESS_raw_data\\ti-2-exsitu\\"
fname1 = fDir + '12\\ff\\ff1_000098.h5'
fname2 = fDir + '12\\ff\\ff2_000098.h5'
yamlFile = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\dex-refined-1.yml"
center = [4888/2,7300/2]
roi_size = 40

# %% Predetermine Cartesian x,y, ,new center, Ainterp
interpDirFile = folder_path + "_interp\\" + folder_path+'_interp_frame_'
all_interp_params = []
for i in range(42000,len(Xs)):
    print(i)
    tth = tths[i]
    eta = etas[i]
    frame = ome_idxs[i]
    Ainterp, new_center, x_cart, y_cart =\
        sf.getInterpParams(center,tth,eta,roi_size,yamlFile)
    interp_params = {"A":Ainterp,"center":new_center,"x_cart":x_cart,"y_cart":y_cart}
    
    try:
        # Load ROI
        roi = sf.loadDexPolarRoi(fname1,fname2,yamlFile,frame,center,tth,eta,roi_size,interp_params)    
        # plt.figure(1)
        # plt.imshow(roi)
        
        # Estimate peak parameters
        tth_vals, eta_vals = np.indices(roi.shape)
        p0 = fitpeak.estimate_pk_parms_2d(eta_vals,tth_vals,roi,"gaussian")
        
        # Fit peak
        p = fitpeak.fit_pk_parms_2d(p0,eta_vals,tth_vals,roi,"gaussian")
        
        # Update tth and eta
        detectDist, mmPerPixel, ff_trans = sf.loadYamlData(yamlFile,center,tth,eta)
        rad_dom, eta_dom = sf.polarDomain(detectDist, mmPerPixel, tth, eta, roi_size)
    
        etaNew = eta_dom[int(np.round(p[1]))]
        radNew = rad_dom[int(np.round(p[2]))]
        tthNew = np.arctan(radNew*mmPerPixel/detectDist)
        
        
        # Update recenter ROI
        Ainterp, new_center, x_cart, y_cart =\
            sf.getInterpParams(center,tthNew,etaNew,roi_size,yamlFile)
        interp_params = {"A":Ainterp,"center":new_center,"x_cart":x_cart,"y_cart":y_cart}
    except:
        print('Recentering spot '+ str(i) + ' failed')
    # roi = sf.loadDexPolarRoi(fname1,fname2,yamlFile,frame,center,tth,eta,roi_size,interp_params)    
    # plt.figure(2)
    # plt.imshow(roi)
    

    all_interp_params.append(interp_params)
    # Save if processed final spot of frame
    if i+1 == len(Xs):
        np.savez( interpDirFile + str(frame) + '.npz',all_interp_params=all_interp_params)
    elif ome_idxs[i+1] > frame:
        np.savez(interpDirFile + str(frame) + '.npz',all_interp_params=all_interp_params)
        all_interp_params = []
# %% Loading data and converting directly to polar coordinates
# spot_ind = 0 
# tth = tths[spot_ind]
# eta = etas[spot_ind]
# frame = ome_idxs[spot_ind]

# b = sf.load_dex(fname1,fname2,frame)
# plt.figure(1)
# plt.imshow(b)
# plt.clim(290,550)
# sf.add_spot_rect(detectDist,mmPerPixel,center,tth,eta,roi_size)


# interpDirFile = folder_path + "interp\\" + folder_path+'_interp_frame_'
# for spot_ind in [0]:
#     # 1. Load spot information
#     tth = tths[spot_ind]
#     eta = etas[spot_ind]
#     frame = ome_idxs[spot_ind]
    
#     # 2. Load interp matrix file if new frame
#     if not frame==ome_idxs[spot_ind]:
#         interp_data = np.load(interpDirFile + str(frame) + '.npz')
#         interp_params = interp_data['all_interp_params']
#         frame = ome_idxs[spot_ind]
    
#     sub_ind = spot_ind - min(np.where(ome_idxs==frame)[0])
#     interp_params[sub_ind]
    
#     # 3. Load ROI and interpolate to polar coordinates 
#     roi_polar = sf.loadDexPolarRoi(fname1,fname2,frame,detectDist,\
#                  mmPerPixel,center,tth,eta,roi_size,ff_trans,interp_params)
        
#     # 4. Compute FWHM with hexrd
#     y_vals, x_vals = np.indices(roi_polar.shape)
#     p0 = fitpeak.estimate_pk_parms_2d(x_vals, y_vals, img,'gaussian')
#     p = fitpeak.fit_pk_parms_2d(p0, x_vals, y_vals, img, 'gaussian')
#     # Add 1e-4 to zero sigma values
#     if p[3]==0: p[3] += 1e-4
#     if p[4]==0: p[4] += 1e-4
#     fit_params_hexrd[i] = p[0:5]
   
#     recon = pkfuncs.gaussian2d(p, x_vals, y_vals)
    
    
# plt.figure(8)
# plt.imshow(roi_polar)
# plt.clim(0,1)