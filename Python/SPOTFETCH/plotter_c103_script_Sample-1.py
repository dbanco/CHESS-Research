#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:14:46 2024

@author: dbanco
"""
import os
import numpy as np
import track_stats_multiplotter as plotter
from multiprocessing import Process
import spotfetch as sf

#### Params for ti-2-tensions data ####
# Output data path
topPath = "/mnt/scratch/dbanco/c103_processing/Sample-1/"
read_path = os.path.join(topPath,'outputs')
spotData = np.load(os.path.join(topPath,'spots','spots.npz'))

# tthdeg = 4
# tth_ring = tthdeg*np.pi/180
# dtth_ring= 2*np.pi/180
# eta1 = np.pi/2
# eta2 = 0
# grains1 = [5]
# spotInds1 = sf.findSpots(spotData,eta=np.pi/2,deta=10,tth=tth_ring,dtth=tth_ring)
# spotInds2 = sf.findSpots(spotData,eta=np.pi,deta=10,tth=tth_ring,dtth=dtth_ring)
# spotInds1 = np.arange(20)
# titleStr1 = f'{eta1} degress'

# grains2 = [11]
# spotInds2 = sf.findSpots(spotData,grains=grains,tth=tth_ring,dtth=dtth_ring)
# spotInds2 = np.arange(20,40)
# titleStr2 = f'{eta2} degress'

grains = [276,288,342]
titleStr = f'Grain {grains[0]}'
spotInds = sf.findSpots(spotData,grains=grains)

processes = []
p1 = Process(target=plotter.start_gui, args=(read_path,spotInds,'Mean',titleStr,spotData,grains))
p1.start()
processes.append(p1)
# p2 = Process(target=plotter.start_gui, args=(read_path, spotInds2, 'Mean',titleStr2,spotData))
# p2.start()
# processes.append(p2)

for p in processes:
   p.join()

