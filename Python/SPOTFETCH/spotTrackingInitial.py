# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 12:16:32 2024

@author: dpqb1
"""

import numpy as np
import spotfetch as sf
from hexrd.fitting import fitpeak
import matplotlib.pyplot as plt

# %% First load indexed spot data
folder_path = "spots_11032023"  
spot_data = np.load(folder_path+'.npz')    
Xs = spot_data['Xs']
Ys = spot_data['Ys']
id_nums = spot_data['id_nums']
tths = spot_data['tths']
etas = spot_data['etas']    
omes = spot_data['omes']    
ome_idxs = spot_data['ome_idxs']  
grain_nums = spot_data['grain_nums']

T = 73
num_spots = len(spot_data['etas'])
tthSeq = np.zeros((num_spots,T+1))
etaSeq = np.zeros((num_spots,T+1))
frmSeq = np.zeros((num_spots,T+1)) 

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
interp_data = np.load(interpDirFile + str(frame) + '.npz', allow_pickle=True)
interp_params = interp_data['all_interp_params']

tthSeq = np.zeros((num_spots,T+1))
etaSeq = np.zeros((num_spots,T+1))
frmSeq = np.zeros((num_spots,T+1)) 
errSeq = np.zeros((num_spots,T+1)) 
pSeq = np.zeros((8,num_spots,T+1)) 

tthSeq[:,0] = tths[spotInds]
etaSeq[:,0] = etas[spotInds]
frmSeq[:,0] = ome_idxs[spotInds]

for t in range(1,T):
    print(t)
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
    
    # for i,ind in enumerate(spotInds[0:1]):
    i = 0
    ind = 0
    tth = tthSeq[i,t-1]
    eta = etaSeq[i,t-1]
    frm = int(frmSeq[i,t-1])

    # Load ROI
    roi = sf.loadDexPolarRoi(fname1,fname2,yamlFile,center,tth,eta,frm,roi_size,interp_params[ind])
    
    plt.figure(1)
    plt.imshow(roi)
    plt.pause(1)
    
    # Estimate peak parameters
    tth_vals, eta_vals = np.indices(roi.shape)
    if np.sum(pSeq[:,i,t-1]) == 0:
        p0 = fitpeak.estimate_pk_parms_2d(eta_vals,tth_vals,roi,"gaussian")
    else:
        p0 = pSeq[:,i,t-1]
    
    # Fit peak
    p = fitpeak.fit_pk_parms_2d(p0,eta_vals,tth_vals,roi,"gaussian")
    residual = fitpeak.fit_pk_obj_2d(p,eta_vals,tth_vals,roi,"gaussian")
    rel_error = np.linalg.norm(residual)/np.linalg.norm(roi.flatten())
    
    detectDist, mmPerPixel, ff_trans = sf.loadYamlData(yamlFile,center,tth,eta)
    rad_dom, eta_dom = sf.polarDomain(detectDist, mmPerPixel, tth, eta, roi_size)
    
    # Spot found if error low enough
    if rel_error < 0.50:
        pSeq[:,i,t] = p
        
        etaSeq[i,t] = eta_dom[int(np.round(p[1]))]
        radNew = rad_dom[int(np.round(p[2]))]
        tthSeq[i,t] = np.arctan(radNew*mmPerPixel/detectDist)
        
        frmSeq[i,t] = frm
        
        errSeq[i,t] = rel_error
    else:
        # Check adjacent omegas for the spot
        adj_frms = [frm-1,frm+1,frm-2,frm+2]
        j = 0
        while rel_error >= 0.50:
            # Load ROI
            roi = sf.loadDexPolarRoi(fname1,fname2,yamlFile,center,tth,eta,adj_frms[j],roi_size,interp_params[ind])
            
            plt.figure(2)
            plt.imshow(roi)
            plt.pause(1)
            # Estimate peak parameters
            tth_vals, eta_vals = np.indices(roi.shape)
            p0 = fitpeak.estimate_pk_parms_2d(eta_vals,tth_vals,roi,"gaussian")
            
            # Fit peak
            p = fitpeak.fit_pk_parms_2d(p0,eta_vals,tth_vals,roi,"gaussian")
            residual = fitpeak.fit_pk_obj_2d(p,eta_vals,tth_vals,roi,"gaussian")
            rel_error = np.linalg.norm(residual)/np.linalg.norm(roi.flatten())
            j += 1
            if rel_error < 0.5:
                pSeq[:,i,t] = p
                
                etaSeq[i,t] = eta_dom[int(np.round(p[1]))]
                radNew = rad_dom[int(np.round(p[2]))]
                tthSeq[i,t] = np.arctan(radNew*mmPerPixel/detectDist)
                
                frmSeq[i,t] = frm
                
                errSeq[i,t] = rel_error
                break
# %% 
plt.figure(1)
plt.plot(etaSeq[0,:])
plt.figure(2)
plt.plot(tthSeq[0,:])
# plt.plot(fwhmEtaSeq[0,:])