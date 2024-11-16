# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:19:11 2024

@author: dpqb1
"""
import numpy as np
import spotfetch as sf
import os
import matplotlib.pyplot as plt

topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
exsituPath = os.path.join(topPath,r"state_0\simulation\outputs\c103_polycrystal_sample_2_state_0_layer_1_output_data.npz")
dataFile = os.path.join(topPath,r"state_*\simulation\outputs\c103_polycrystal_sample_2_state_*_layer_1_output_data.npz")

# 1. Spot data npz files
spotFiles = []
fullscanRange = np.arange(0,5)
for state in fullscanRange:
    spotsOut = os.path.join(topPath,f'state_{state}')
    spotFiles.append(os.path.join(spotsOut,'spots.npz'))

# 2. Load data
spotData = []
for i in range(len(spotFiles)):
    spotData.append(np.load(spotFiles[i]))


# 3. Get number of spots with ome_idxs == 0
for i in range(len(spotFiles)):
    num_spots = sum(spotData[i]['grain_nums'] == 0)
    display(num_spots)

num_spots = 150
eta_diffs = np.zeros((num_spots,4))
tth_diffs = np.zeros((num_spots,4))

grain_num = 0
for j in range(num_spots):
    for i in range(len(spotFiles)-1):     
        sel_grain = spotData[i]['grain_nums'] == grain_num
        eta_ji = spotData[i]['etas'][sel_grain][j]
        tth_ji = spotData[i]['tths'][sel_grain][j]
        
        sel_grain = spotData[i+1]['grain_nums'] == grain_num
        eta_jip1 = spotData[i+1]['etas'][sel_grain][j]
        tth_jip1 = spotData[i+1]['tths'][sel_grain][j]
        
        eta_diffs[j,i] = abs(eta_jip1 - eta_ji)
        tth_diffs[j,i] = abs(tth_jip1 - tth_ji)
   
plt.figure
for i in range(4):
    plt.subplot(4,1,i+1)
    plt.hist(eta_diffs[:,i])
    
# for i in range(len(spotFiles)):
#     num_spots = sum(spotData[i]['id_nums'] == 3)
#     display(num_spots)
    
    
# for i in range(len(spotFiles)):
#     num_spots = sum(spotData[i]['ome_idxs'] == 2)
#     display(num_spots)

# 3. Compute tth, eta differences in time 
# etaDiffs = []
# tthDiffs = []
# for i in range(len(spotFiles)-1):
#     etaDiffs.append(spotData[i+1]['etas'] - spotData[i]['etas'])
#     tthDiffs.append(spotData[i+1]['tths'] - spotData[i]['tths']) 
    
# i = 0
# spotData[i+4]['etas'][0] 
# spotData[i+3]['etas'][0] 
# spotData[i+2]['etas'][0] 
# spotData[i+1]['etas'][0] 
# spotData[i]['etas'][0]

# plt.hist(etaDiffs[0])

