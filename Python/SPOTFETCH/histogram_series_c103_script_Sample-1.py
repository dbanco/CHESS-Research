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
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pickle
from hexrd.fitting import peakfunctions as pkfuncs

#### Params for ti-2-tensions data ####
# Output data path
topPath = "C:\\Users\\dpqb1\\Documents\\Data\\c103_2024"
# topPath = "/nfs/chess//user/dbanco/c103_processing/Sample-1/"
read_path = os.path.join(topPath, 'outputs')

spotsDir = "C:\\Users\\dpqb1\\Documents\\Data\\c103_processing"
spotData = np.load(os.path.join(spotsDir, 'spots.npz'))
# spotData = np.load(os.path.join(topPath,'spots','spots.npz'))

grains = [44, 158]
tth110 = 0.06562438
tth200 = 0.09267698
tth211 = 0.1136209
tth321 = 0.1736603
tthhkl = [tth110, tth211]
hklNames = ['110', '211']

spotIndsList = []
titleStrList = []
for grain in grains:
    for i, tth in enumerate(tthhkl):
        spotIndsList.append(sf.findSpots(
            spotData, grains=[grain], tth=tth, dtth=0.01))
        titleStrList.append(f'Grain {grain}, {hklNames[i]} : $FWHM_\eta$')
sInds = sf.findSpots(spotData, grains=grains)

# Plot a series of histograms for each grain
trackData = []
for i, spotInds in enumerate(spotIndsList):
    trackData.append([])
    for k in spotInds:
        with open(os.path.join(read_path, f'trackData_{k}.pkl'), 'rb') as f:
            trackData[i].append(pickle.load(f))

scanRange = np.concatenate(
    (np.array([368, 372, 376, 380]), np.arange(383, 406), [407]))
T = len(scanRange)
K1 = len(spotIndsList)

dataCounts = []
dataBins = []
numPlots = len(spotIndsList)
dataArray = []

# Update all features
for ii, spotInds in enumerate(spotIndsList):
    dataArray.append(np.zeros((6, len(spotInds), T)))
    dataArray[ii][:] = np.nan
    for k in range(len(spotInds)):
        if k > len(trackData[ii])-1:
            continue
        for t in range(T):
            if t > len(trackData[ii][k])-1:
                continue
            if len(trackData[ii][k][t]) > 0:

                avgFWHMeta, avgFWHMtth, avgEta, avgTth = sf.compAvgParams(
                    trackData[ii][k][t], 0.4)
                FWHMome = sf.estFWHMomega(trackData[ii][k][t])
                Ome = sf.estMEANomega(trackData[ii][k][t])
                dataArray[ii][0, k, t] = FWHMome
                dataArray[ii][2, k, t] = avgFWHMeta
                dataArray[ii][4, k, t] = avgFWHMtth
                dataArray[ii][1, k, t] = Ome
                dataArray[ii][3, k, t] = avgEta
                dataArray[ii][5, k, t] = avgTth

fig, axes = plt.subplots(nrows=3, ncols=int(np.ceil(T/3)), figsize=(T, 16))

histTime = np.zeros((12, T, numPlots))
for ii in range(numPlots):
    for t in range(T):
        ax = axes.ravel()[t]
        hist_t = ax.hist(dataArray[ii][2, :, t], bins=12, range=(0, 1.25))
        histTime[:, t, ii] = hist_t[0]
        ax.set_title(f'Histogram at Time Step {t+1}')
        ax.set_xlabel('Value')
        ax.set_ylabel('Frequency')

bin_centers = 0.5 * (hist_t[1][1:] + hist_t[1][:-1])


plt.tight_layout()
plt.show()

vmin = histTime.min()
vmax = histTime.max()

fig2, axes2 = plt.subplots(nrows=2, ncols=2, figsize=(30, 10))
fSize = 14
for ii in range(numPlots):
    im1 = axes2.ravel()[ii].imshow(histTime[:, :, ii])
    cbar2 = fig2.colorbar(im1, ax=axes2.ravel()[ii])
    cbar2.ax.tick_params(labelsize=fSize)
    axes2.ravel()[ii].set_title(titleStrList[ii], fontsize=fSize)
    axes2.ravel()[ii].set_yticks(np.arange(len(bin_centers)))
    axes2.ravel()[ii].set_yticklabels(
        [f'{center:.2f}' for center in bin_centers], fontsize=fSize)
    axes2.ravel()[ii].set_xlabel('Time index', fontsize=fSize)
    axes2.ravel()[ii].set_ylabel('Bin centers $(^\circ)$', fontsize=fSize)

plt.tight_layout()
plt.show()

# Plot each ROI, Recon, and Recon err
fig3, axes3 = plt.subplots(nrows=10, ncols=20, figsize=(50, 50))
j = 0
for ii, spotInds in enumerate(spotIndsList):
    for k in range(len(spotInds)):
        if k > len(trackData[ii])-1:
            continue
        for t in range(T):
            if t > len(trackData[ii][k])-1:
                continue
            if len(trackData[ii][k][t]) > 0:
                for i in range(len(trackData[ii][k][t])):
                    roi = trackData[ii][k][t][i]['roi']
                    tth_vals, eta_vals = np.indices(roi.shape)
                    f = pkfuncs.gaussian2d(
                        trackData[ii][k][t][i]['p'], eta_vals, tth_vals)
                    im1 = axes3.ravel()[j].imshow(roi)
                    axes3.ravel()[j].set_xticks([])
                    axes3.ravel()[j].set_yticks([])
                    j += 1
                    im2 = axes3.ravel()[j].imshow(f)
                    errf = trackData[ii][k][t][i]['err']
                    axes3.ravel()[j].set_title(f'{errf:1.3}')
                    axes3.ravel()[j].set_xticks([])
                    axes3.ravel()[j].set_yticks([])
                    j += 1
                    if j == 200:
                        break
            if j == 200:
                break
        if j == 200:
            break
    if j == 200:
        break
        # if j == 200:
        #     fig3, axes3 = plt.subplots(nrows=10, ncols=20, figsize=(50, 50))
        #     j = 0

spotErrs = []
for ii, spotInds in enumerate(spotIndsList):
    for k in range(len(spotInds)):
        if k > len(trackData[ii])-1:
            continue
        for t in range(T):
            if t > len(trackData[ii][k])-1:
                continue
            if len(trackData[ii][k][t]) > 0:
                for i in range(len(trackData[ii][k][t])):
                    err_t = trackData[ii][k][t][i]['err']
                    if err_t < 0.4:
                        spotErrs.append(err_t)

# Distribution of spot errors
plt.figure()
plt.hist(np.array(spotErrs), bins=20, range=(0, 1))

len(spotErrs)
# spotErrs[spotErrs<0.5] =

# Histogram to view distribution of spots by radius (ring)
plt.figure()
tths = spotData['tths']
plt.hist(tths, bins=50, range=(0, 0.17))


# import h5py
# hf = h5py.File("C:\\Users\\dpqb1\\Documents\\Data\\c103_processing\\c103_sam2_refined.h5", 'r')
# data = hf.dataset

# import chess_detectors as cd
# yf = cd.read_yaml("C:\\Users\\dpqb1\\Documents\\Data\\c103_processing\\c103_processing\\eiger16M_monolith_mruby_062224_FINAL.yml")
# yf['detectors']


# plotter.start_gui(read_path, spotIndsList, 'Mean',titleStr,spotData,grains)

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
