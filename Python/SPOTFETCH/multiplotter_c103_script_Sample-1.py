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
topPath = "C:\\Users\\dpqb1\\Documents\\Data\\c103_2024"
# topPath = "/nfs/chess//user/dbanco/c103_processing/Sample-1/"
read_path = os.path.join(topPath,'outputs')

spotsDir = "C:\\Users\\dpqb1\\Documents\\Data\\c103_processing"
spotData = np.load(os.path.join(spotsDir,'spots.npz'))
#spotData = np.load(os.path.join(topPath,'spots','spots.npz'))

grains = [44,158]
tth110 = 0.06562438
tth200 = 0.09267698
tth211 = 0.1136209
tth321 = 0.1736603
tthhkl = [tth110,tth200,tth211]

spotIndsList = []
for grain in grains:
    # for i,tth in enumerate(tthhkl):
    spotIndsList.append(sf.findSpots(spotData,grains=[grain],tth=tth211,dtth=0.01))

titleStr = f'Grain {grains}'

plotter.start_gui(read_path, spotIndsList, 'Mean',titleStr,spotData,grains)

# dome = 3
# scanRange = np.concatenate((np.array([364,368,372,376,380]), np.arange(383,406), [407]))
# trackPath = os.path.join(topPath,'outputs')
# sf.roiTrackVisual(spotIndsList[0],spotData,dome,scanRange,trackPath,dataPath,params):

# processes = []
# p1 = Process(target=plotter.start_gui, args=(read_path, spotIndsList, 'Mean',titleStr,spotData,grains))
# p1.start()
# processes.append(p1)

# for p in processes:
#    p.join()

