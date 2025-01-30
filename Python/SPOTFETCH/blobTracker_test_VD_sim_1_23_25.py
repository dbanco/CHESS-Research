# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 08:53:32 2024

@author: dpqb1
"""
import numpy as np
import spotfetch as sf
import os

# topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
topPath = r"E:\VD_sim_processing"
dataFile = os.path.join(topPath,r"state_{scan}\simulation\outputs\c103_polycrystal_sample_2_state_{scan}_layer_1_output_data.npz")

scanRange = np.arange(5)
spotsFiles = []
for state in scanRange:
    spotsDir = os.path.join(topPath,f'state_{state}')
    spotsFiles.append(os.path.join(spotsDir,"spots.npz"))

prmBlob = {}
prmBlob['detector'] = 'eiger_sim'
prmBlob['imSize'] = (5000,5000)
prmBlob['yamlFile'] = os.path.join(topPath,'c103_eiger_calibration.yml') #mruby_0401_eiger_calibration
prmBlob['roiSize'] = [31,31,11]
prmBlob['initDist'] = 4
prmBlob['gamma'] = [8,8,8,6,2,2]
prmBlob['pool'] = 16
prmBlob['parallelFlag'] = False
prmBlob['benchmarkFlag'] = True

prmBlob['gaussian_derivatives'] = False
prmBlob['sigmas'] = np.array([2,2,2])
prmBlob['dsigma'] = np.array([1.5,1.5,1.5])

spotInd = 0
suffix = '_1_23_25_blob_'
ttPath = os.path.join(topPath,'outputs'+suffix)
if not os.path.exists(ttPath):
    os.mkdir(ttPath)

dataFileSequence = [] 
for scan in scanRange:
    dataFileSequence.append(dataFile.format(scan=scan)) 
    
initSpotData = np.load(spotsFiles[0])

spotInds = np.arange(1)

# Run tracking
for spotInd in spotInds:
    sf.trackSpotBlob(spotInd,initSpotData,dataFileSequence,ttPath,prmBlob)

# Make track image files
output_path = os.path.join(topPath,'imageFigs'+suffix)
dome = 5
num_cols = 5
sf.makeBlobTrackImages(dome,num_cols,output_path,spotInds,initSpotData,scanRange,dataFileSequence,ttPath,spotsFiles,prmBlob)
