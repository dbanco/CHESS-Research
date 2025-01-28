# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 13:34:37 2025

Script to load in .pkl file containing a single 
roiTensor with shape [M_tth,M_eta,M_ome,T]

The ROI is extracted at the same [tth,eta] location at each omega and each scan.
In the data I have provided, M_ome = 7 and the spot of interest should at the 
very least appear centered in frame roiTensor[:,:,3,0] 

HEXD spot data, shared here has no spot tracking information at all. It is possible
for the spot to exit the ROI either in eta or omega. Spot 4650 for example, the
first spot, does not appear to move.


@author: Daniel Banco
"""
import pickle
import os 
import matplotlib.pyplot as plt
import numpy as np

# Collapse this function for plotting
def plotROI(roiTensor,num_figs):
    T = roiTensor.shape[3]
    M_ome = roiTensor.shape[2]
    
    num_cols = int(np.ceil(T/num_figs))
    scanRange = np.arange(29)
    
    fig_list = []
    axes_list = []
    
    # Create a figure and a set of subplots
    for i in range(num_figs):
        fig, axes = plt.subplots(M_ome, num_cols, figsize=(20, 15))
        
        # Remove x and y ticks for clarity
        for ax_row in axes:
            for ax in ax_row:
                ax.set_xticks([])
                ax.set_yticks([])
            
        fig_list.append(fig)
        axes_list.append(axes)
        
        i1 = 0 + i*num_cols
        i2 = num_cols-1 + i*num_cols
        j1 = scanRange[i1]
        if i2 >= T:
            j2 = scanRange[-1]
        else:
            j2 = scanRange[i2]
        # Add common labels
        fig_list[i].text(0.04, 0.5, r'$\omega$ frame', va='center', rotation='vertical', fontsize=24)
        fig_list[i].text(0.5, 0.04, 'Scan #', ha='center', fontsize=24)
        fig_list[i].text(0.5, 0.95, f'Spot {spotInd}, Scans {j1}-{j2}', ha='center', fontsize=32)      
        fig_list[i].subplots_adjust(wspace=0.05, hspace=0.01)
        
    for scan_ind in range(T):
        for om_ind in range(M_ome):
            i = int(np.floor(scan_ind/num_cols))
            j = np.mod(scan_ind,num_cols)
            ax = axes_list[i][om_ind,j]
            scan = scanRange[scan_ind]
            
            # 1. Show ROI
            roi = roiTensor[:,:,om_ind,scan_ind]
            ax.imshow(roi)
            ax.text(1, 2, f'max: {roi.max():.2f}', color='white', fontsize=12, weight='bold')
          
            # Label plots
            if om_ind == M_ome-1:
                ax.set_xlabel(f'{scan}')
        

# Load and plot ROI
spotInd = 4650
topDir = r'E:\Data\c103_processing\roiTensors_grain_44'
roiFile = os.path.join(topDir,f'roiTensor_{spotInd}.pkl')
with open(roiFile, 'rb') as f:
    roiTensor = pickle.load(f)
    
num_figs = 2 # Number figures across which to show roi data
plotROI(roiTensor,num_figs)

