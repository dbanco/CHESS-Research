# -*- coding: utf-8 -*-
"""
Created on Thu May  2 09:04:54 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
from hexrd.fitting import fitpeak
import matplotlib.pyplot as plt

# %% First load indexed spot data
folder_path = "spots_11032023"  
grain_data = np.load(folder_path+'.npz')    
Xs = grain_data['Xs']
Ys = grain_data['Ys']
id_nums = grain_data['id_nums']
tths = grain_data['tths']
etas = grain_data['etas']    
omes = grain_data['omes']    
ome_idxs = grain_data['ome_idxs']  
grain_nums = grain_data['grain_nums']

T = 73
tthSeq = np.zeros((len(Xs),T+1))
etaSeq = np.zeros((len(Xs),T+1))
frmSeq = np.zeros((len(Xs),T+1)) 

# %% intial time step Dataset Parameters
fDir = "D:\\CHESS_raw_data\\ti-2-exsitu\\"
fname1 = fDir + '12\\ff\\ff1_000098.h5'
fname2 = fDir + '12\\ff\\ff2_000098.h5'
yamlFile = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\dex-refined-1.yml"


# %% Process Spot Data Tracking in Omega over time
interpDirFile = folder_path + "_interp\\" + folder_path+'_interp_frame_'
numSpots = 46181
center = [4888/2,7300/2]
roi_size = 40

# Process only spots that begin on frame 4 to start
frame = 4
spotInds = np.where(ome_idxs == frame)[0]

ampSeq = np.zeros((len(spotInds),T+1))
tthSeq = np.zeros((len(spotInds),T+1))
etaSeq = np.zeros((len(spotInds),T+1))
frmSeq = np.zeros((len(spotInds),T+1)) 
fwhmEtaSeq = np.zeros((len(spotInds),T+1)) 
fwhmTthSeq = np.zeros((len(spotInds),T+1)) 
errSeq = np.zeros((len(spotInds),T+1)) 

tthSeq[:,0] = tths[spotInds]
etaSeq[:,0] = etas[spotInds]
frmSeq[:,0] = ome_idxs[spotInds]

t = 1

dirNum = t+1
# Do stupid logic because we are missing a timestep
if dirNum == 62: 
    t +=1
    dirNum += 1
if dirNum < 62:
    fileNum = 105 + dirNum
else:
    fileNum = 105 + dirNum - 1

dNum = str(dirNum)
fNum = str(fileNum)
fDir = "D:\\CHESS_raw_data\\ti-2-tension\\"
fname1 = fDir + dNum + '\\ff\\ff1_000' + fNum + '.h5'
fname2 = fDir + dNum + '\\ff\\ff2_000' + fNum + '.h5'

interp_data = np.load(interpDirFile + str(frame) + '.npz', allow_pickle=True)
interp_params = interp_data['all_interp_params']

# for i,ind in enumerate(spotInds[0:1]):
i = 0
ind = 0
tth = tthSeq[i,t-1]
eta = etaSeq[i,t-1]
frm = int(frmSeq[i,t-1])

# Load ROI
roi = sf.loadDexPolarRoi(fname1,fname2,yamlFile,frm,center,tth,eta,roi_size,interp_params[ind])

plt.figure(1)
plt.imshow(roi)
plt.pause(1)

# Estimate peak parameters
tth_vals, eta_vals = np.indices(roi.shape)
p0 = fitpeak.estimate_pk_parms_2d(eta_vals,tth_vals,roi,"gaussian")

# Fit peak
p = fitpeak.fit_pk_parms_2d(p0,eta_vals,tth_vals,roi,"gaussian")
residual = fitpeak.fit_pk_obj_2d(p,eta_vals,tth_vals,roi,"gaussian")
rel_error = np.linalg.norm(residual)/np.linalg.norm(roi.flatten())

# Spot found if error low enough



