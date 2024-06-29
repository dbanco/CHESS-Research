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
topPath = "/nfs/chess//user/dbanco/c103_processing/Sample-1/"
read_path = os.path.join(topPath,'outputs')

spotData = np.load(os.path.join(topPath,'spots','spots.npz'))

grains = [276,288,342]
spotIndsList = []
for grain in grains:
    spotIndsList.append(sf.findSpots(spotData,grains=[grain]))

# spotInds1 = sf.findSpots(spotData,grains=grains[0],tth=tth_ring,dtth=dtth_ring)
# spotInds2 = sf.findSpots(spotData,grains=grains[1],tth=tth_ring,dtth=dtth_ring)
# spotInds1 = np.arange(20)
# spotInds2 = np.arange(20,40)
# spotIndsList = [spotInds1,spotInds2]

titleStr = f'Grain {grains}'

processes = []
p1 = Process(target=plotter.start_gui, args=(read_path, spotIndsList, 'Mean',titleStr,spotData,grains))
p1.start()
processes.append(p1)

for p in processes:
   p.join()

