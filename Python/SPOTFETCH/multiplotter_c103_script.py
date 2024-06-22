#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:14:46 2024

@author: dbanco
"""
import os
import numpy as np
import spotfetch as sf
import track_stats_multiplotter as plotter
from multiprocessing import Process

#### Params for ti-2-tensions data ####
# Output data path
topPath = "/mnt/scratch/dbanco/c103_processing/Sample-1/layer-1"
read_path = os.path.join(topPath,'outputs')

spotsDir = "spots_11032023"
spotsFile = spotsDir + ".npz"  
spotData = np.load(os.path.join(topPath,'spots',spotsFile))

tthdeg = 3.76
tth_ring = tthdeg*np.pi/180
dtth_ring= 4*np.pi/180

grains = [5,11]
spotIndsList = []
for grain in grains:
    spotIndsList.append(sf.findSpots(spotData,grains=[grain],\
                                              tth=tth_ring,dtth=dtth_ring))

# spotInds1 = sf.findSpots(spotData,grains=grains[0],tth=tth_ring,dtth=dtth_ring)
# spotInds2 = sf.findSpots(spotData,grains=grains[1],tth=tth_ring,dtth=dtth_ring)
# spotInds1 = np.arange(20)
# spotInds2 = np.arange(20,40)
# spotIndsList = [spotInds1,spotInds2]

titleStr = f'Grain {grains} at {tthdeg} degress'

processes = []
p1 = Process(target=plotter.start_gui, args=(read_path, spotIndsList, 'Mean',titleStr,spotData,grains))
p1.start()
processes.append(p1)


for p in processes:
   p.join()

