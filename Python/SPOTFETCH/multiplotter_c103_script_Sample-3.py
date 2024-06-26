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
topPath = "/mnt/scratch/dbanco/c103_processing/Sample-3/"
read_path = os.path.join(topPath,'outputs')
spotData = np.load(os.path.join(topPath,'spots','spots.npz'))

grains = [1,2,3]

tthdeg = 4
tth_ring = tthdeg*np.pi/180
dtth_ring= 2*np.pi/180
deta = 20*np.pi/180

spotInds1 = sf.findSpots(spotData,eta=0,deta=deta,tth=tth_ring,dtth=tth_ring)
spotInds2 = sf.findSpots(spotData,eta=np.pi/4,deta=deta,tth=tth_ring,dtth=dtth_ring)
spotInds3 = sf.findSpots(spotData,eta=np.pi/2,deta=deta,tth=tth_ring,dtth=dtth_ring)

spotIndsList = [spotInds1,spotInds2,spotInds3]

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

