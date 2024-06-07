# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:13:58 2024

@author: dpqb1
"""

import numpy as np
import matplotlib.pyplot as plt
# import os
import spotfetch as sf
from hexrd.fitting import fitpeak
from hexrd.fitting import peakfunctions as pkfuncs
import time
# import h5py


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
detectDist = 608
mmPerPixel = 0.0748
center = [4888/2,7300/2]
ff1_tx = 111.77123331882476
ff1_ty = -0.8903409731826509
ff2_tx = -129.61824519824643
ff2_ty = -1.9743054988931037
ff_trans = [ff1_tx,ff1_ty,ff2_tx,ff2_ty]

roi_size = 40

# %% Testing efficient loading and FWHM processing of polar ROIs
# spot_ind = 0 
# tth = tths[spot_ind]
# eta = etas[spot_ind]
# frame = ome_idxs[spot_ind]

# b = sf.load_dex(fname1,fname2,frame)
# plt.figure(1)
# plt.imshow(b)
# plt.clim(290,550)
# sf.add_spot_rect(detectDist,mmPerPixel,center,tth,eta,roi_size)

# %%
interpDirFile = folder_path + "_interp\\" + folder_path+'_interp_frame_'
FWHM_etas = np.zeros(etas.shape) 
numSpots = 46181
computeTimes = np.zeros(etas.shape) 
old_frame = -1
for spot_ind in range(numSpots):
    print(spot_ind)
    tic = time.time()
    # 1. Load spot information
    tth = tths[spot_ind]
    eta = etas[spot_ind]
    frame = ome_idxs[spot_ind]
    # 2. Load interp matrix file if new frame
    if (not frame == old_frame) or (not 'interp_params' in locals()):
        interp_data = np.load(interpDirFile + str(frame) + '.npz', allow_pickle=True)
        interp_params = interp_data['all_interp_params']     
    
    sub_ind = spot_ind - min(np.where(ome_idxs==frame)[0])
    
    toc1 = time.time()
    
    # 3. Load ROI and interpolate to polar coordinates 
    roi_polar = sf.loadDexPolarRoi(fname1,fname2,frame,detectDist,\
                 mmPerPixel,center,tth,eta,roi_size,ff_trans,interp_params[sub_ind])
        
    toc2 = time.time()
        
    # 4. Compute FWHM with hexrd
    try:
        y_vals, x_vals = np.indices(roi_polar.shape)
        p0 = fitpeak.estimate_pk_parms_2d(x_vals, y_vals, roi_polar,'gaussian')
        p = fitpeak.fit_pk_parms_2d(p0, x_vals, y_vals, roi_polar, 'gaussian')
        # Add 1e-4 to zero sigma values
        if p[3]==0: p[3] += 1e-4
        if p[4]==0: p[4] += 1e-4
        FWHM_etas[spot_ind] = p[3]   
        toc = time.time()
        computeTimes[spot_ind] = toc-tic
    except:
        print('Frame was empty?Size 0?')
    old_frame = frame
    # recon = pkfuncs.gaussian2d(p, x_vals, y_vals)
    outStr = 'Time1: ' + str(toc1-tic) +\
        ', Time2: ' + str(toc2-toc1) +\
        ', Time3: ' + str(toc-toc2) +\
        ', Total: ' + str(toc-tic)
    print(outStr)