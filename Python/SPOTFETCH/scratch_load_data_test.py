# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
import os
import time
import pickle

# %% Processing Setup

# Output data path
topPath_nfs = "/nfs/chess/user/dbanco/c103_processing/Sample-1/layer-1"
topPath_mnt = "/mnt/scratch/dbanco/c103_processing/Sample-1/layer-1"

# Raw HEXD Data
exsituPath = "/nfs/chess/raw/2024-2/id3a/miller-3528-c/c103-1-ff-1/1/ff"
dataPath = "/mnt/scratch/pagan-dwellfatigue-series-cycle-2024-2/"
initTracksPath_nfs = os.path.join(topPath_nfs,'initTracks')
initTracksPath_mnt = os.path.join(topPath_mnt,'initTracks')


# Full dexela image size and roi size
params = {}

params['detector'] = 'eiger'
params['imSize'] = (5000,5000) 
params['yamlFile'] = '/mnt/scratch/dbanco/mruby_eiger_calibration_single_grain_v01.yml'

params['roiSize'] = 40
params['gamma'] = [8,8,4,4]

newData = False
scan = 10
while True:
    # Try reading in file for new scan
    fname = os.path.join(dataPath,f'ti6242-sample1-1_00{scan}_EIG16M_CdTe_000001.h5')
    try:
        img = sf.load_eiger(fname,params,0)
        img1 = img.copy()
        sum1d = np.sum(img,0)
        outFile = os.path.join(topPath_mnt,'outputs',f'sum1d_{scan}.pkl')
        with open(outFile, 'wb') as f:
            pickle.dump(sum1d, f)
            print('Wrote data')
            
    except:
        print('File not read')
        time.sleep(1)
        scan = 10
    
    
    scan += 1