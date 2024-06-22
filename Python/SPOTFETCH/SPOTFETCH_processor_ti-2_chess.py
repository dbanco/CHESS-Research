# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:29 2024

@author: dpqb1
"""
# %% Processing Setup #####
import numpy as np
import spotfetch as sf
import os

#   1. Set up paths: exsitu, data, outputs
exsituPath = "/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-exsitu/12/ff/"
topPath = "/nfs/chess/user/dbanco/ti-2_processing"
dataPath = "/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/"

# %% 2. Load in or collect spots data
spotsDir = "spots_11032023"
# Get data from spots data
# sf.collectSpotsData(topPath, spotsDir)
spotsFile = spotsDir + ".npz"  
spotData = np.load(os.path.join(topPath,'spots',spotsFile))

# %% 3. Detector and tracking parameters
params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = os.path.join(topPath,'dex-refined-1.yml')
params['imSize'] = (4888,7300) 
params['roiSize'] = 40
params['gamma'] = [3,5,4,4] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16

# %% 4. Inspect spot tracks on initial scan data
scanNum = 1
fNum = 106
frame = 10
fname1 = os.path.join(dataPath,f'{scanNum}','ff','ff1_000' + str(fNum) + '.h5')
fname2 = os.path.join(dataPath,f'{scanNum}','ff','ff2_000' + str(fNum) + '.h5')
sf.plotSpotWedges(spotData,fname1,fname2,frame,params)

# %% 5. Determine which spots to process
# spotInds = np.arange(120,1000)
grains = [5,11]
tth_ring = 3.76*np.pi/180
dtth_ring= 4*np.pi/180
spotInds = sf.findSpots(spotData,grains=grains,tth=tth_ring,dtth=dtth_ring)

# %% 6. Begin Processing
scan1 = 43
initTracksPath = os.path.join(topPath,'outputs')
sf.initExsituTracks(initTracksPath,exsituPath, spotData, spotInds, params, scan1-1)
sf.spotTracker(dataPath,topPath,exsituPath,spotData,spotInds,params,scan1)
