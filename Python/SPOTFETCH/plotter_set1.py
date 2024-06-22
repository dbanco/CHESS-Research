#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:14:46 2024

@author: dbanco
"""
import os
import threading
import numpy as np
import spotfetch as sf
import SPOTFETCH_multiplotter as plotter

#### Params for ti-2-tensions data ####
# Output data path
topPath = "/nfs/chess/user/dbanco/ti-2_processing"
spotsDir = "spots_11032023"
spotsFile = spotsDir + ".npz"  
spotData = np.load(os.path.join(topPath,'spots',spotsFile))

read_path1 = "/nfs/chess/user/dbanco/ti-2_processing/outputs"
grains = [5]
tth_ring = 3.76*np.pi/180
dtth_ring= 4*np.pi/180
spotInds = sf.findSpots(spotData,grains=grains,tth=tth_ring,dtth=dtth_ring)
titleStr = 'Grain 5 at 3.76 degress'
plotter.start_gui(read_path1, spotInds, 'Mean',titleStr)


'''
#### Params for c103 Sample-1 ####
topPath = "/mnt/scratch/...."
spotsFile = os.path.join(topPath,'spots','spots.npz')
spotData = np.load(spotsFile)
read_path = os.path.join(topPath,'outputs')

# Set 1
spotInds = sf.findSpots(spotData,5,np.pi/2,0.1)
plotter.start_gui(read_path,spotInds,'Delta')

# Set 2
spotInds = sf.findSpots(spotData,50,np.pi/2,0.1)
plotter.start_gui(read_path,spotInds,'Delta')
'''
