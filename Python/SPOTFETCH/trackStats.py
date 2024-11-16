# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 10:28:24 2024

@author: dpqb1
"""
import numpy as np
import os
import spotfetch as sf
import pickle
import matplotlib.pyplot as plt

#### Params for ti-2-tensions data ####
# Output data path
topPath = "C:\\Users\\dpqb1\\Documents\\Data\\c103_processing\\c103_processing\\Sample-1"
# topPath = "/nfs/chess//user/dbanco/c103_processing/Sample-1/"
read_path = os.path.join(topPath, 'outputs')

spotData = np.load(os.path.join(topPath, 'spots', 'spots.npz'))

grains = [276, 288, 342]
spotIndsList = []
for grain in grains:
    spotIndsList.append(sf.findSpots(spotData, grains=[grain]))

plt.figure()
# Load in track data
spotTracks = []
for i, spotInds in enumerate(spotIndsList):
    print(i)
    spotTracks.append(np.zeros(len(spotInds)))
    for j, k in enumerate(spotInds):
        file_path = os.path.join(read_path, f'trackData_{k}.pkl')
        with open(file_path, 'rb') as f:
            trackData = pickle.load(f)
            trackCount = 0
            for track in trackData:
                omLen = len(track)
                if omLen > 0:
                    trackCount += 1
            spotTracks[i][j] = trackCount

    plt.subplot(1, 3, i+1)
    plt.hist(spotTracks[i])
    plt.xlabel('# of spots detected (0-30)')
    plt.ylabel('Counts')
    plt.title(f'Grain {grains[i]}')
