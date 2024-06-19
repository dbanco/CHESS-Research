# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 10:19:33 2024

@author: dpqb1
"""
import spotfetch as sf
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os

topPath = "/nfs/chess/user/dbanco/ti-2_processing"
outPath = os.path.join(topPath,'outputs')
dataPath = "/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/"

# Get initial tracks
spotsFile= os.path.join(topPath,'spots','spots_11032023.npz')
spotData = np.load(spotsFile)

initData = {
    'tths': spotData['tths'],
    'etas': spotData['etas'],
    'frms': spotData['ome_idxs']
}

# Choose spot
for k in [8]:#range(10):
    print(f'Showing Spot {k}')
    eta0 = initData['etas'][k]
    tth0 = initData['tths'][k]
    frm0 = initData['frms'][k]
    
    # Define Omega range
    om0 = frm0
    omMin = om0-4
    omMax = om0+4
    omRange = np.arange(omMin,omMax+1)
    for i,om in enumerate(omRange):
        omRange[i] = sf.wrapFrame(om)
    
    # Time range
    scanRange = np.arange(43,62)
    
    # Define the number of rows and columns
    rows = len(omRange)
    columns = len(scanRange)
    
    # Create a figure and a set of subplots
    fig, axes = plt.subplots(rows, columns, figsize=(20, 15))
    
    # Remove x and y ticks for clarity
    for ax_row in axes:
        for ax in ax_row:
            ax.set_xticks([])
            ax.set_yticks([])
    
    # Add common labels
    fig.text(0.04, 0.5, r'$\omega$ frame', va='center', rotation='vertical', fontsize=24)
    fig.text(0.5, 0.04, 'Scan #', ha='center', fontsize=24)
    fig.text(0.5, 0.95, f'Spot {k}', ha='center', fontsize=32)
    
    # Path to the track data file
    track_file = os.path.join(outPath,f'trackData_{k}.pkl')
    
    # Full dexela image size and roi size
    params = {}
    # params['detector'] = 'eiger'
    params['detector'] = 'dexela'
    params['yamlFile'] = "/nfs/chess/user/dbanco/ti-2_processing/dex-refined-1.yml"
    params['imSize'] = (4888,7300) 
    params['roiSize'] = 40
    
    # Function to update plots
    def update_plots(track_data):
        
        # Organize all tracks
        T = len(track_data)
        omTrack = np.zeros((rows,columns)) - 1
        scanTrack = np.zeros((rows,columns)) - 1
        etaRoiTrack = np.zeros((rows,columns)) + eta0
        tthRoiTrack = np.zeros((rows,columns)) + tth0
        etaTrack = np.zeros((rows,columns)) + eta0
        tthTrack = np.zeros((rows,columns)) + tth0
        p1Track = np.zeros((rows,columns))
        p2Track = np.zeros((rows,columns))
        
        for t in range(T):
            track = track_data[t]
            if len(track) > 0:
                for j in range(len(track)):
                    scan = track[j]['scan']
                    frm = track[j]['frm']
                    
                    if (frm in omRange) & (scan in scanRange):
                        ind1 = int(np.where(scanRange == scan)[0][0])
                        ind2 = int(np.where(omRange == frm)[0][0])
                    
                        eta = track[j]['eta']
                        tth = track[j]['tth']
                        etaRoi = track[j]['etaRoi']
                        tthRoi = track[j]['tthRoi']
                        p = track[j]['p']
                        
                        scanTrack[ind2,ind1] = scan
                        omTrack[ind2,ind1] = frm
                        etaTrack[ind2,ind1] = eta
                        tthTrack[ind2,ind1] = tth
                        etaRoiTrack[ind2,ind1] = etaRoi
                        tthRoiTrack[ind2,ind1] = tthRoi
                        p1Track[ind2,ind1] = p[1]
                        p2Track[ind2,ind1] = p[2]
                    
        for j in range(columns):
            fname1, fname2 = sf.timeToFile(scanRange[j],dataPath)            
            for i in range(rows): 
                ax = axes[i, j]
                ax.clear()
    
                etaRoi = etaRoiTrack[i,j]
                tthRoi = tthRoiTrack[i,j]
                frm = omRange[i]
                scan = scanRange[j]
                p1 = p1Track[i,j]
                p2 = p2Track[i,j]
                # Show roi
                roi = sf.loadPolarROI(fname1,fname2,tthRoi,etaRoi,frm,params)
                ax.imshow(roi)
                
                if (p1 > 0) & (p2 > 0):
                    # Plot rect
                    eta = etaTrack[i,j]
                    tth = tthTrack[i,j]
                    detectDist, mmPerPixel, ff_trans = sf.loadYamlData(params,tth,eta)
                    boxSize = 10
                    start_col = round(p1 - boxSize//2)
                    start_row = round(p2 - boxSize//2)
                    rect = plt.Rectangle((start_col, start_row), boxSize, boxSize,
                                          linewidth=1, edgecolor='r', facecolor='none')
                    ax.add_patch(rect)
                    
                if i == rows-1:
                    ax.set_xlabel(f'{scan}')
                if j == 0:
                    ax.set_ylabel(f'{frm}')
                    
                ax.set_xticks([])
                ax.set_yticks([])
        plt.draw()

    # Initial plot setup
    if os.path.exists(track_file):
        with open(track_file, 'rb') as f:
            track_data = pickle.load(f)
            update_plots(track_data)
