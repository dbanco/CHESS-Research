# -*- coding: utf-8 -*-
"""
Created on Thu May  2 10:07:45 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
# from hexrd.fitting import fitpeak
# import matplotlib.pyplot as plt

# %% First load indexed spot data
folder_path = "spots_11032023"  
spot_data = np.load(folder_path+'.npz')    
yamlFile = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\dex-refined-1.yml"
interpDirFile = folder_path + "_interp\\" + folder_path+'_interp_frame_'

numSpots = len(spot_data['etas'])
imSize = (4888,7300)
roi_size = 40


# %% Plot range of spots
# Notes:
#   ind = 0, frame =4, timeFrms = range(53,68), omeFrms = range(2,7)

ind = 0
tth = spot_data['tths'][ind]
eta = spot_data['etas'][ind]
frame = 4
timeFrms = range(51,61)
omeFrms = range(-2,9)
# timeFrms = range(43,45)
# omeFrms = range(3,6)

sf.roiAdjacent(yamlFile,ind,tth,eta,imSize,roi_size,frame,omeFrms,timeFrms)

            