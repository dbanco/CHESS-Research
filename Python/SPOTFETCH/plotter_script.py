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
spotInds1 = sf.findSpots(spotData,5,np.pi/2,0.1)
# spotInds2 = np.arange(120,1000)
spotInds2 = np.arange(8)
threading.Thread(target=plotter.start_gui, args=(read_path1, spotInds2, 'Delta'), daemon=True).start()
threading.Thread(target=plotter.start_gui, args=(read_path1, spotInds2, 'Mean'), daemon=True).start()

# spotInds3 = np.arange(8,16)
# threading.Thread(target=plotter.start_gui, args=(read_path1, spotInds3, 'Delta'), daemon=True).start()
# threading.Thread(target=plotter.start_gui, args=(read_path1, spotInds3, 'Mean'), daemon=True).start()
 

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
