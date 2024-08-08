# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 14:06:56 2024

@author: dpqb1
"""
import matplotlib.pyplot as plt
import spotfetch as sf
import numpy as np
import pickle
import time

topPath = "C:\\Users\\dpqb1\\Documents\\Data\\ti-2\\"
outPath = "C:\\Users\\dpqb1\\Documents\\Data\\ti-2\\outputs\\"
dataPath = "D:\\CHESS_raw_data\\ti-2-tension\\"
# Define the number of rows and columns
rows = 7
columns = 7
k = 0

# Create a figure and a set of subplots
fig, axes = plt.subplots(rows, columns, figsize=(14, 10))  # Adjust figsize as needed

# Get initial tracks
spotsDir = "spots_11032023"
spotsFile = spotsDir + ".npz" 
spotData = np.load(topPath + spotsFile)
spotInds = np.arange(0,4)

# Full dexela image size and roi size
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = "C:\\Users\\dpqb1\\Documents\\Data\\ti-2\\dex-refined-1.yml"
params['imSize'] = (4888,7300) 
params['roiSize'] = 40

initData = {}
initData['tths'] = spotData['tths'][spotInds]
initData['etas'] = spotData['etas'][spotInds]
initData['frms'] = spotData['ome_idxs'][spotInds]

# Load in current track for spot k
outFile = outPath + f'trackData_{k}.pkl'
with open(outFile, 'rb') as f:
    trackData = pickle.load(f)
    
om0 = trackData[0][0]['frm']
omMin = om0-3
omMax =om0+3

scan0 = trackData[0][0]['scan']

eta = initData['etas'][spotInds[k]]
tth = initData['tths'][spotInds[k]]
frm = initData['frms'][spotInds[k]]

scn = scan0 + 11
times = np.arange(scn,scn+7)
oms = np.arange(omMin,omMax+1)
for i,t in enumerate(times):
    fname1, fname2 = sf.timeToFile(t,dataPath)
    track = trackData[i]
    for j,om in enumerate(oms):
        # Load Roi at Omega j
        frm =  sf.wrapFrame(om)
        roi = sf.loadPolarROI(fname1,fname2,tth,eta,frm,params)
            
        # Plot example data on each subplot
        axes[j, i].imshow(roi)
        
        # Add rectangle for track
        
        # Remove x and y ticks for clarity (optional)
        axes[i, j].set_xticks([])
        axes[i, j].set_yticks([])
    for j in range(len(track)):
        frm = track[j]['frm']
        tth = track[j]['tth']
        eta = track[j]['eta']
        jj = np.where(oms == frm)[0][0]
        
        # 1. Load YAML data
        detectDist, mmPerPixel, ff_trans = sf.loadYamlData(params,tth,eta)
        
        # 2. Add rectangle
        boxSize = 10
        x = track[j]['p'][1]
        y = track[j]['p'][2]
        start_col = round(x - boxSize//2)
        start_row = round(y - boxSize//2)
        rect = plt.Rectangle((start_col, start_row), boxSize, boxSize,
                              linewidth=1, edgecolor='r', facecolor='none')
        axes[jj,i].add_patch(rect)
        
# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()


# Save the figure (optional)
# plt.savefig("subplots_5x7.png", dpi=300)

# Show the figure
plt.show()
