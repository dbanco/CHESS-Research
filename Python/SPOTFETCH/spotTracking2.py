# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:27:15 2024

@author: dpqb1
"""
import numpy as np
import spotfetch as sf
import matplotlib.pyplot as plt
import pickle
import os

# %% Processing Setup

# Output data path
output_path = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\outputs"

# Initial indexing data
folder_path = "spots_11032023"  
spot_data = np.load(folder_path+'.npz')

# HEXD Data
fDir = "D:\\CHESS_raw_data\\ti-2-tension\\"

# Detector parameters
yamlFile = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\dex-refined-1.yml"

# Initial interpolation matrix
# interpDirFile = folder_path + "_interp\\" + folder_path+'_interp_frame_'

# Initialize Spotfetch tracking parameters
tInds = np.arange(43,62)
spotInds = np.arange(0,2)
# frame = 4
# spotInds = np.where(spot_data['ome_idxs'] == frame)[0]

params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\dex-refined-1.yml"
params['imSize'] = (4888,7300)
params['roiSize'] = 40

trackData, initData = sf.trackData(spot_data,spotInds)
currTracks = []
for i,t in enumerate(tInds):
    print('')
    print(f'Time: {t}')
    fname1,fname2 = sf.timeToFile(t,fDir)
    trackData.append([])
    print('Spot:', end=" ")
    for j,s in enumerate(spotInds):#initData['etas'].shape[0]):    
        print(f'{j}', end=" ")
        trackData[i].append([])
        if i == 0:           
            eta = initData['etas'][s]
            tth = initData['tths'][s]
            frm = initData['frms'][s]
            prevTracks = []
            newTrack, peakFound = sf.evaluateROI(fname1,fname2,prevTracks,\
                                tth,eta,int(frm),t,params)
            trackData[i][j].append(newTrack)
        else:
            prevTracks = trackData[i-1][j]
            if len(prevTracks) == 0:
                continue
            
        # Initial Search: through all current omega tracks, then check up and down for\
        # tracks (not sure exactly when search through omega will be considered done)    
        for track in prevTracks:
            eta = track['eta']
            tth = track['tth']
            frm = track['frm']
            # print(f'Checking prev track at frame {frm}')
            # Load ROI and fit peak
            newTrack, peakFound = sf.evaluateROI(fname1,fname2,prevTracks,\
                                tth,eta,int(frm),t,params)
            
            # Add to list if peakFound
            if peakFound: 
                # print(f'Peak found at frame {frm}')
                trackData[i][j].append(newTrack)
                 
        # Next search up and down in omega for more peaks detections and append 
        # to front or back of list depending on where a track is found. If there
        # are no tracks to begin with, we should know where the last known track
        # was to be able to search omegas we have not checked yet. There will
        # have to be an arbitrary limit on the omegas to be checked
        #
        # If we have tracks, we can just stop searching in omega once the first
        # fails to be detected in either direction
        #
        # The criterion for detecting a peak should be looser when we have no
        # current track and it can be tighter when we actually have a track
        
        # Get last known time index for track
        # tIdx = i
        # while len(trackData[s,tIdx]) == 0:
        #     tIdx -= 1
        
        # If we have a track
        if len(trackData[i][j]) > 0:
            compareTrack = trackData[i][j]
            frm1 = trackData[i][j][0]['frm']
            frm2 = trackData[i][j][-1]['frm']
        else:
            compareTrack = trackData[i-1][j]
            frm1 = trackData[i-1][j][0]['frm']
            frm2 = trackData[i-1][j][-1]['frm']
        
        # Search down
        count = 0
        while count < 3:
            frm1 = sf.wrapFrame(frm1 - 1)
            
            # Load ROI and fit peak
            newTrack, peakFound = sf.evaluateROI(fname1,fname2,compareTrack,\
                                tth,eta,int(frm1),t,params)
  
            # Add to list if peakFound
            if peakFound: 
                print(f'Found more at {frm1}')
                trackData[i][j].insert(0,newTrack)
                count = 0
            else:
                count += 1
 
        # Search Up
        count = 0
        while count < 3:
            frm2 = sf.wrapFrame(frm2 + 1)
            # Load ROI and fit peak
            newTrack, peakFound = sf.evaluateROI(fname1,fname2,compareTrack,\
                                tth,eta,int(frm2),t,params)
                       
            # Add to list if peakFound
            if peakFound: 
                # print(f'Found more at {frm2}')
                trackData[i][j].append(newTrack)
                count = 0
            else:
                count += 1
    
    """Write processed data to a JSON file named trackData.json."""
    with open(os.path.join(output_path,'trackData.pkl'), 'wb') as f:
        pickle.dump(trackData, f)



# sf.scatterOmeg(trackData)

# # %%
# # sf.plotTrackData(trackData,spotInds)
# num_plots = 1#trackData['pSeq'].shape[1]
# num_rows = 6
# num_cols = 5

# fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_rows, num_cols))

# # Plot each image in a subplot
# for i in range(num_plots):
#     ax = axes.flat[i]
#     ax.plot(tInds,trackData['etaSeq'][i,:])
#     ax.set_title(f'Spot {i} $\mu_\eta$')
#     # ax.set_xlabel('Time')
    
    
# fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_rows, num_cols))

# # Plot each image in a subplot
# for i in range(num_plots):
#     ax = axes.flat[i]
#     ax.plot(tInds,trackData['pSeq'][3,i,:])
#     ax.set_title(f'Spot {i} $FWHM_\eta$')
#     # ax.set_xlabel('Time')
    
# fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_rows, num_cols))

# # Plot each image in a subplot
# for i in range(num_plots):
#     ax = axes.flat[i]
#     ax.plot(tInds,trackData['pSeq'][0,i,:])
#     ax.set_title(f'Spot {i} $Amplitude$')
#     # ax.set_xlabel('Time')
    
# fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_rows, num_cols))

# # Plot each image in a subplot
# for i in range(num_plots):
#     ax = axes.flat[i]
#     ax.plot(tInds,trackData['frmSeq'][i,:])
#     ax.set_title(f'Spot {i} $OmegFrmSeq$')
#     # ax.set_xlabel('Time')
    
# fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_rows, num_cols))

# # Plot each image in a subplot
# for i in range(num_plots):
#     ax = axes.flat[i]
#     ax.plot(tInds,trackData['errSeq'][i,:])
#     ax.set_title(f'Spot {i} $Error$')
#     # ax.set_xlabel('Time')

# # Need code to automatically fullscreen and save each one of these

