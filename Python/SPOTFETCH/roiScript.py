# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 22:11:18 2024

@author: dpqb1
"""
import numpy as np 
import os
import yaml
from matplotlib import pyplot as plt
from hexrd import xrdutil
from hexrd import transforms


topPath = r'/nfs/chess/user/seg246/software/MechVD/VD_simulations/c103_polycrystal_sample_2'
spotsFile = os.path.join(topPath,'state_0','spots.npz')
spotData = np.load(spotsFile)
fname = os.path.join(topPath,r"state_0\simulation\outputs\c103_polycrystal_sample_2_state_0_layer_1_output_data.npz")

# Load eiger calib yaml
yamlFile = os.path.join(topPath,'c103_eiger_calibration.yml')
with open(yamlFile, 'r') as file:
    yamlData = yaml.safe_load(file)
trans = yamlData['detectors']['eiger']['transform']['translation']
tilt = yamlData['detectors']['eiger']['transform']['tilt']
detectDist = -trans[2]
mmPerPixel = yamlData['detectors']['eiger']['pixels']['size'][0]

# Load image frame_i (2-1441)
frame_i = 2
simData = np.load(fname)       
shp = simData['shape']
rowD = simData[f'{frame_i-2}_row']
colD = simData[f'{frame_i-2}_col']
datD = simData[f'{frame_i-2}_data']
img = np.zeros((shp[0],shp[1]))
for i in range(len(rowD)):
    img[rowD[i],colD[i]] = datD[i]

# Show image
fig = plt.figure()
# b[b>100] = 100
plt.imshow(img)
plt.clim(np.median(img),np.median(img)+20)
plt.colorbar()

# Inputs to project spot locationsonto detector plane
rMat_d = transforms.xf.makeDetectorRotMat(tilt)
rMat_c = np.eye(3)#transforms.xfcapi.makeRotMatOfExpMap()
chi = 0
tVec_d = np.array(trans)
tVec_c = np.zeros((3,1))
tVec_s = np.zeros((3,1))
distortion = None #np.eye(3)

# Get indices of spots at frame_i
frms = spotData['ome_idxs']
spotInds = np.where(frms == frame_i)[0]

# Use translation/tilt to plot spot locations on image at frame_i
for k in spotInds:
    # Spot location
    tth = spotData['tths'][k]
    eta = spotData['etas'][k]
    ome = spotData['omes'][k]
    allAngs = np.array([[tth,eta,ome]])
    # Project on detector
    det_xy, rmats_s, on_plane = xrdutil._project_on_detector_plane(allAngs,
                                                                rMat_d, rMat_c, chi,
                                                                tVec_d, tVec_c, tVec_s,
                                                                distortion)
    # Convert from real detector space to pixel space
    det_pix = det_xy.ravel()/mmPerPixel
    frame_i = spotData['ome_idxs'][k]
    row = det_pix[0] + (shp[1]-1)/2
    col = -det_pix[1] + (shp[0]-1)/2
    # Plot over image
    plt.plot(row,col,'x',markersize=12)

