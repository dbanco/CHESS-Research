#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:14:46 2024

@author: dbanco
"""
import os
import numpy as np
import track_stats_plotter as plotter
from multiprocessing import Process

#### Params for ti-2-tensions data ####
# Output data path
topPath = "/mnt/scratch/dbanco/c103_processing/Sample-1/layer-1"
read_path = os.path.join(topPath,'outputs')

spotData = np.load(os.path.join(topPath,'spots','spots.npz'))

tthdeg = 3.76
tth_ring = tthdeg*np.pi/180
dtth_ring= 1*np.pi/180

grains1 = [5]
# spotInds1 = sf.findSpots(spotData,grains=grains1,tth=tth_ring,dtth=dtth_ring)
spotInds1 = np.arange(20)
titleStr1 = f'Grain {grains1[0]} at {tthdeg} degress'

grains2 = [11]
# spotInds2 = sf.findSpots(spotData,grains=grains,tth=tth_ring,dtth=dtth_ring)
spotInds2 = np.arange(20,40)
titleStr2 = f'Grain {grains2[0]}at {tthdeg} degress'

processes = []
p1 = Process(target=plotter.start_gui, args=(read_path, spotInds1, 'Mean',titleStr1,spotData))
p1.start()
processes.append(p1)
p2 = Process(target=plotter.start_gui, args=(read_path, spotInds2, 'Mean',titleStr2,spotData))
p2.start()
processes.append(p2)

for p in processes:
   p.join()

