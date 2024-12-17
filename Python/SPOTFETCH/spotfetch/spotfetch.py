# -*- coding: utf-8 -*-
"""
spotfetch.py

Main module for the spotfetch package: tools for X-ray diffraction data analysis 
and spot tracking in imaging experiments.

Created on: Fri Apr 12 13:22:47 2024
Author: Daniel Banco
Email: dpqb10@gmail.com
Version: 0.1.0

Description:
This module serves as the main entry point for the spotfetch package, providing 
tracking and analysis tools for X-ray diffraction data. 
It integrates submodules for data processing, detector interfacing, data 
labeling and visualization.

License:
MIT License

Copyright (c) 2024 dpqb1

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import numpy as np
import warnings
import pandas as pd
import os
import matplotlib.pyplot as plt
from hexrd.fitting import fitpeak
import glob
from .detectors import loadPolarROI,loadYamlData,polarDomain

FRAME1 = 2
NUMFRAMES = 1440
OMEG_RANGE = 360

def omegaToFrame(omega,startFrame=FRAME1,numFrames=NUMFRAMES,omegaRange=OMEG_RANGE,startOmeg = 0):
    step = omegaRange/numFrames
    frame = (np.floor((omega*180/np.pi-startOmeg)/step) + startFrame).astype(int)
    return frame
      
def frameToOmega(frame,startFrame=FRAME1,numFrames=NUMFRAMES,omegaRange=OMEG_RANGE,startOmeg = 0):
    step = omegaRange/numFrames
    omega = (frame - startFrame)*step + startOmeg
    return omega

def mapOmega(omega):
    if omega > 180:
        omega = omega - 360
    return omega

def mapDiff(diff):
    return (diff + 180) % 360 - 180

def wrapFrame(frm,frm0=FRAME1,numFrms=NUMFRAMES):
    return np.mod(frm-frm0,numFrms)+frm0

def pathToFile(path):
    fnames = os.listdir(path)
    fnames.sort()
    for i in range(len(fnames)):
        fnames[i] = os.path.join(path,fnames[i])
    return fnames

def timeToFile(t,fDir):
    dNum = str(t)
    topDir = os.path.join(fDir + dNum, 'ff')
    fnames = pathToFile(topDir)
    return fnames

def getDataFileSequence(dataFile,scanRange):
    dataFileSequence = []
    for scan in scanRange:
        template = dataFile.format(num2=scan)
        pattern = os.path.join(template)
        fname = glob.glob(pattern)[0]
        if len(fname) > 0:
            dataFileSequence.append(fname)
    return dataFileSequence

def findSpots(spotData, **kwargs):
    grains = kwargs.get('grains', None)
    tths = kwargs.get('tths', None)
    dtth = kwargs.get('dtth', None)
    eta = kwargs.get('eta', None)
    deta = kwargs.get('deta', None)
    frm = kwargs.get('frm', None)

    
    if grains != None:
        grain_nums = spotData['grain_nums']
        cond1 = np.array(grain_nums.shape)
        cond1 = False
        for g in grains:
            cond1 = cond1 | (grain_nums == g)
    else:
        cond1 = True
    if eta != None:
        etas = spotData['etas']
        cond2 = (etas > eta-deta) & (etas < eta+deta)
    else:
        cond2 = True
    if tths != None:
        tths = spotData['tths']
        cond3 = np.array(tths.shape)
        cond3 = False
        for tth in tths:
            tths = spotData['tths']
            cond3 = cond3 | (tths > tth-dtth) & (tths < tth+dtth)
    else:
        cond3 = True
    if frm != None:
        frms = spotData['ome_idxs']
        cond4 = frms == frm
    else:
        cond4 = True
    
    # if cond1 & cond2 & cond3 & cond4:
    #     grain_nums = spotData['grain_nums']
    #     spotInds = np.ones(grain_nums.shape)
    #     spotInds = np.where(spotInds)[0]
    # else:
    spotInds = np.where(cond1 & cond2 & cond3 & cond4)[0]
        
    return spotInds

def collectSpotsData(outPath,spotsPath,sortFlag=False):
    # Get a list of file names in the directory and sort them
    all_entries = os.listdir(spotsPath)
    directories = [entry for entry in all_entries if os.path.isdir(os.path.join(spotsPath,entry))]
    fold_names = sorted(directories)
    
    created = 0
    # Loop through sorted file names
    for fold_name in fold_names:
        file_names = sorted(os.listdir(os.path.join(spotsPath, fold_name)))
        print(fold_name)
        for file_name in file_names:
            if file_name.endswith(".out"):  # Check if the file has a ".out" extension
                file_path = os.path.join(spotsPath,fold_name, file_name)
                
                # Load .out file
                df = pd.read_csv(file_path,sep=' \s+',engine='python')  
                
                # Get a HEXD spot location
                if not created:
                    Xs = np.array(df['pred X'])
                    Ys = np.array(df['pred Y'])
                    Xm = np.array(df['meas X'])
                    Ym = np.array(df['meas Y'])
                    id_nums = np.array(df['# ID'])
                    PID = np.array(df['PID'])
                    H = np.array(df['H'])
                    K = np.array(df['K'])
                    L = np.array(df['L'])
                    grain_nums = int(file_name[-7:-4])*np.ones(df['pred X'].shape)
                    tths = np.array(df['meas tth'])
                    etas = np.array(df['meas eta'])
                    omes = np.array(df['meas ome'])
                    tths_pred = np.array(df['pred tth'])
                    etas_pred = np.array(df['pred eta'])
                    omes_pred = np.array(df['pred ome'])
                    created = 1
                else:
                    Xs = np.append(Xs, np.array(df['pred X']))
                    Ys = np.append(Ys, np.array(df['pred Y']))
                    Xm = np.append(Xm, np.array(df['meas X']))
                    Ym = np.append(Ym, np.array(df['meas Y']))                    
                    id_nums = np.append(id_nums, np.array(df['# ID']))
                    PID = np.append(PID, np.array(df['PID']))
                    H = np.append(H, np.array(df['H']))
                    K = np.append(K, np.array(df['K']))
                    L = np.append(L, np.array(df['L']))
                    grain_nums = np.append(grain_nums,int(file_name[-7:-4])*np.ones(df['pred X'].shape))
                    tths = np.append(tths, np.array(df['meas tth']))
                    etas = np.append(etas, np.array(df['meas eta']))
                    omes = np.append(omes, np.array(df['meas ome']))
                    tths_pred = np.append(tths_pred, np.array(df['pred tth']))
                    etas_pred = np.append(etas_pred, np.array(df['pred eta']))
                    omes_pred = np.append(omes_pred, np.array(df['pred ome']))
                
    invalid_nums = (id_nums == -999) | np.isnan(tths.astype(float)) | np.isnan(etas.astype(float))       
    Xs = np.delete(Xs, invalid_nums).astype(float)
    Ys = np.delete(Ys, invalid_nums).astype(float)
    Xm = np.delete(Xm, invalid_nums).astype(float)
    Ym = np.delete(Ym, invalid_nums).astype(float)
    id_nums = np.delete(id_nums, invalid_nums).astype(float)
    PID = np.delete(PID, invalid_nums).astype(float)
    H = np.delete(H, invalid_nums).astype(float)
    K = np.delete(K, invalid_nums).astype(float)
    L = np.delete(L, invalid_nums).astype(float)
    grain_nums = np.delete(grain_nums, invalid_nums).astype(float)
    tths = np.delete(tths, invalid_nums).astype(float)
    etas = np.delete(etas, invalid_nums).astype(float)
    omes = np.delete(omes, invalid_nums).astype(float)
    tths_pred = np.delete(tths_pred, invalid_nums).astype(float)
    etas_pred = np.delete(etas_pred, invalid_nums).astype(float)
    omes_pred = np.delete(omes_pred, invalid_nums).astype(float)
    
    ome_idxs = omegaToFrame(omes)
    
    if sortFlag:
        sort_ind = np.argsort(omes)
        ome_idxs = ome_idxs[sort_ind]
        Xs = Xs[sort_ind]
        Ys = Ys[sort_ind]
        Xm = Xm[sort_ind]
        Ym = Ym[sort_ind]
        id_nums = id_nums[sort_ind]
        PID = PID[sort_ind]
        H = H[sort_ind]
        K = K[sort_ind]
        L = L[sort_ind]
        tths = tths[sort_ind]
        etas = etas[sort_ind]
        omes = omes[sort_ind]
        tths_pred = tths_pred[sort_ind]
        etas_pred = etas_pred[sort_ind]
        omes_pred = omes_pred[sort_ind]
        grain_nums = grain_nums[sort_ind]
    saveFile = os.path.join(outPath,'spots.npz')
    np.savez(saveFile,Xs=Xs,Ys=Ys,Xm=Xm,Ym=Ym,id_nums=id_nums,
    tths=tths,etas=etas,omes=omes,
    tths_pred=tths_pred,etas_pred=etas_pred,omes_pred=omes_pred,
    ome_idxs=ome_idxs,grain_nums=grain_nums,PID=PID,H=H,K=K,L=L)

def getSpotID(spotData,k):
    grNum = spotData['grain_nums'][k]
    PID = spotData['PID'][k]
    H = spotData['H'][k]
    K = spotData['K'][k]
    L = spotData['L'][k]
    idNum = spotData['id_nums'][k]
    spot_id = [grNum,PID,H,K,L,idNum]
    return spot_id

def matchSpotID(spotData,spot_id,k):
    grNum = spotData['grain_nums']
    PID = spotData['PID']
    H = spotData['H']
    K = spotData['K']
    L = spotData['L']
    idNum = spotData['id_nums']
    
    bin_array = (grNum == spot_id[0]) & (PID == spot_id[1]) & \
                (H == spot_id[2]) & (K == spot_id[3]) & (L == spot_id[4])
    match_ind = np.where(bin_array)[0]
    
    if len(match_ind) > 1:
        diffs = abs(idNum[match_ind] - spot_id[5])
        match_ind = match_ind[diffs<50]
        diffs = diffs[diffs<50]
        if len(match_ind) > 1:
            ind = np.argmin(diffs)
            match_ind = match_ind[ind]
            
    if (not np.isscalar(match_ind)) and (len(match_ind) == 0):
        print(f'Match not found for spot {k}')
        return np.nan
    
    return match_ind

def estMEANomega(track):
    step = OMEG_RANGE/NUMFRAMES
    if len(track) == 1:
        return frameToOmega(track[0]['frm'])
    
    roiOmega = np.zeros(len(track))
    omegaRange = np.zeros(len(track))
    meanOmega = 0
    for i in range(len(track)):
        roiOmega[i] = np.nansum(track[i]['roi'],axis=None)
        omegaRange[i] = frameToOmega(track[i]['frm'])
        meanOmega += i*roiOmega[i]
    meanOmega = meanOmega/np.sum(roiOmega)
    
    ind1 = int(np.floor(meanOmega))
    ind2 = int(np.ceil(meanOmega))
    
    if omegaRange[ind1] > omegaRange[ind2]:
        meanOmega = omegaRange[ind1] + (meanOmega-ind1)*step
    else:
        meanOmega = meanOmega*step

    # Map from 180-360 to -180-0 if first omega greater than last omega 
    # if track[-1]['frm'] > track[0]['frm']:
    meanOmega = mapOmega(meanOmega)
    return meanOmega

def estFWHMomega(track):
    step = OMEG_RANGE/NUMFRAMES
    if len(track) == 1:
        return 0
    
    roiOmega = np.zeros(len(track))
    omegaRange = np.zeros(len(track))
    meanOmega = 0
    varOmega = 0
    for i in range(len(track)):
        roiOmega[i] = np.nansum(track[i]['roi'],axis=None)
        omegaRange[i] = frameToOmega(track[i]['frm'])
        meanOmega += i*roiOmega[i]
    meanOmega = meanOmega/np.nansum(roiOmega)
    
    for i in range(len(track)):
        varOmega += roiOmega[i]*(i-meanOmega)**2/np.nansum(roiOmega,None)*step**2

    fwhmOmega = 2*np.sqrt(2*np.log(2)*varOmega)
    
    return fwhmOmega
    
def loadSpotsAtFrame(spot_data,fnames,frame,params,detectFrame=[]):
    tths = spot_data['tths']
    etas = spot_data['etas']     
    ome_idxs = spot_data['ome_idxs'] 
    
    # Get spot indices at frame
    spotInds = np.where(ome_idxs == frame)[0]
    
    # Init ROIs list so they can be different sizes
    roi_list = []
    
    for i, ind in enumerate(spotInds):
       
        # 1. Load spot information
        tth = tths[ind]
        eta = etas[ind]

        # 2. Load polar ROI
        if detectFrame == []:
            detectFrame = frame
            
        roi_polar = loadPolarROI([fnames],tth,eta,detectFrame,params)    
        
        roi_list.append(roi_polar)
        
    return roi_list

def fitModel(eta_vals,tth_vals,roi,params):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if params['peak_func'] == "gaussian":
                p0 = fitpeak.estimate_pk_parms_2d(eta_vals,tth_vals,roi,"gaussian")
                p = fitpeak.fit_pk_parms_2d(p0,eta_vals,tth_vals,roi,"gaussian")
                if (p[3]==0): p[3] += 0.001
                if (p[4]==0): p[4] += 0.001
            elif params['peak_func'] == "gaussian_rot":
                p0 = fitpeak.estimate_pk_parms_2d(eta_vals,tth_vals,roi,"gaussian_rot")
                p = fitpeak.fit_pk_parms_2d(p0,eta_vals,tth_vals,roi,"gaussian_rot")
                if (p[3]==0): p[3] += 0.001
                if (p[4]==0): p[4] += 0.001
            peakFound = True
            # Make sure peak lies within ROI
            roiSize = params['roiSize']
            if (p[1] > roiSize[0]-0.5) | (p[2] > roiSize[1]-0.5) | (p[1] < -0.5) | (p[2] < -0.5):
                peakFound = False
            return p, peakFound
    except:
        peakFound = False
        return 0, peakFound 
     
def visualTrackData(trackData):
    T = len(trackData)
    K = len(trackData[0])
    FWHMeta = np.zeros((T,K))
    meaneta = np.zeros((T,K))
    
    for t in range(T):
        if trackData[t] != []:
            for k in range(K):
                if trackData[t][k] != []:
                    FWHMeta[t] = trackData[t][k][0]['p'][3]
                    meaneta[t] = trackData[t][k][0]['eta']
    return FWHMeta,meaneta
            
def compAvgParams(track,errorThresh):   
    J = len(track) 
    avgFWHMeta = 0
    avgFWHMtth = 0
    avgEta = 0
    avgTth = 0
    numJ = 0
    for j in range(J):
        if track[j]['err'] < errorThresh:
            avgFWHMeta += track[j]['p'][3]*track[j]['deta']
            avgFWHMtth += track[j]['p'][4]*track[j]['dtth']
            avgEta += track[j]['eta']
            avgTth += track[j]['tth']
            numJ += 1            
    if numJ > 0:
        avgEta = avgEta*(180/np.pi)/numJ
        avgTth = avgTth*(180/np.pi)/numJ
        avgFWHMeta = avgFWHMeta*(180/np.pi)/numJ
        avgFWHMtth = avgFWHMtth*(180/np.pi)/numJ
    else:
        avgEta = np.nan
        avgTth = np.nan
        avgFWHMeta = np.nan
        avgFWHMtth = np.nan
    return avgFWHMeta,avgFWHMtth,avgEta,avgTth

def showROI(ax,dataPath,scan,frm,tthRoi,etaRoi,params):
    
    roi = loadROI(dataPath,scan,frm,etaRoi,tthRoi,params)
        
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlData(params,tthRoi,etaRoi)
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel,tthRoi, etaRoi, params['roiSize'])
    tth_dom = np.arctan(rad_dom*mmPerPixel/detectDist)
    
    ax.imshow(roi)
    plt.xticks([0,params['roiSize'][0]], [f'{eta_dom[0]:.4g}',f'{eta_dom[-1]:.4g}'])
    plt.yticks([0,params['roiSize'][1]], [f'{tth_dom[0]:.4g}',f'{tth_dom[-1]:.4g}'])
  
    ax.text(1, 2, f'max: {roi.max():.2f}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 4, f'scan: {scan}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 6, f'frame: {frm}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 4, f'eta: {etaRoi:.4f}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 6, f'tth: {tthRoi:.4f}', color='white', fontsize=12, weight='bold')

def showROIpolar(ax,dataPath,scan,frm,tthRoi,etaRoi,params):
    
    roi = loadROI(dataPath,scan,frm,etaRoi,tthRoi,params)
        
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlData(params,tthRoi,etaRoi)
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel,\
                          tthRoi, etaRoi, params['roiSize'])
    tth_dom = np.arctan(rad_dom*mmPerPixel/detectDist)
    
    ax.imshow(roi)
    plt.xticks([0,39], [f'{eta_dom[0]:.4g}',f'{eta_dom[-1]:.4g}'])
    plt.yticks([0,39], [f'{tth_dom[0]:.4g}',f'{tth_dom[-1]:.4g}'])
  
    ax.text(1, 2, f'max: {roi.max():.2f}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 4, f'scan: {scan}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 6, f'frame: {frm}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 4, f'eta: {etaRoi:.4f}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 6, f'tth: {tthRoi:.4f}', color='white', fontsize=12, weight='bold')

def showInitial(ax,etaRoi,tthRoi,eta0,tth0,params):
    # Hexrd location
    [y_pos, x_pos] = etaTthToPix(eta0,tth0,etaRoi,tthRoi,params)
    if (y_pos >= 0) and (y_pos < params['roiSize'][0]) and\
       (x_pos >= 0) and (x_pos < params['roiSize'][1]): 
        
        ax.plot(x_pos,y_pos,marker='s',markersize=10,\
                fillstyle='none',color='white')

def showTrack(ax,track,etaRoi,tthRoi,eta0,tth0,params):
    eta = track['eta']
    tth = track['tth']
    etaPrev = track['etaRoi']
    tthPrev = track['tthRoi']
    # Current track
    [y_pos, x_pos] = etaTthToPix(eta,tth,etaRoi,tthRoi,params)
    if (y_pos >= 0) and (y_pos < params['roiSize'][0]) and\
       (x_pos >= 0) and (x_pos < params['roiSize'][1]):
        # ax.plot(x_pos,y_pos,marker='s',markersize=8,\
        #         fillstyle='none',color='red')
        ax.plot(x_pos,y_pos,marker='o',markersize=13,fillstyle='none',color='red')
    # Previous track
    # [y_pos, x_pos] = etaTthToPix(etaPrev,tthPrev,etaRoi,tthRoi,params)
    # if (y_pos >= 0) and (y_pos < params['roiSize'][0]) and\
    #    (x_pos >= 0) and (x_pos < params['roiSize'][1]):
    #     ax.plot(x_pos,y_pos,marker='o',markersize=10,\
    #             fillstyle='none',color='yellow')


def showTruth(ax,truth,etaRoi,tthRoi,params):
    eta1 = truth['eta1']
    tth1 = truth['tth1']
    eta2 = truth['eta2']
    tth2 = truth['tth2']
    [y1, x1] = etaTthToPix(eta1,tth1,etaRoi,tthRoi,params)
    [y2, x2] = etaTthToPix(eta2,tth2,etaRoi,tthRoi,params)
    
    rect = plt.Rectangle((x1, y1), x2-x1, y2-y1,
                          linewidth=1, edgecolor='g', facecolor='none')
    ax.add_patch(rect)

def xyToEtaTth(x,y,params):
    detectDist, mmPerPixel, trans, tilt = loadYamlData(params)
    
    eta = np.arctan2(y,x)
    rad = np.sqrt(x**2 + y**2)
    tth = np.arctan(rad/detectDist)
    
    return eta, tth

def xyToEtaTthRecenter(x,y,params):
    detectDist, mmPerPixel, trans, tilt = loadYamlData(params)
    
    x = x + trans[0]
    y = y + trans[1]
    
    eta = np.arctan2(y,x)
    rad = np.sqrt(x**2 + y**2)
    tth = np.arctan(rad/detectDist)
    
    return eta, tth
    
def etaTthToPix(eta,tth,etaRoi,tthRoi,params):
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlData(params,tthRoi,etaRoi)
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel,\
                          tthRoi, etaRoi, roiSize)
        
    rad = np.tan(tth)*detectDist/mmPerPixel
    
    row_pos = (rad-rad_dom[0])/(rad_dom[-1]-rad_dom[0])*(roiSize[0]-1)
    col_pos = (eta-eta_dom[0])/(eta_dom[-1]-eta_dom[0])*(roiSize[1]-1)
    
    return row_pos, col_pos

def pixToEtaTth(p1,p2,tthRoi,etaRoi,params):
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlData(params,tthRoi,etaRoi)
    rad_dom, eta_dom = polarDomain(detectDist,mmPerPixel,tthRoi,etaRoi,roiSize)  
    
    deta = abs(eta_dom[1] - eta_dom[0])
    hypot = detectDist*np.cos(tthRoi)
    dtth = np.arctan(mmPerPixel/hypot)
    
    i1 = int(np.floor(p1))
    j1 = int(np.floor(p2))
    
    etaNew = eta_dom[i1] + deta*np.mod(p1,1)
    radNew = rad_dom[j1] + np.mod(p2,1)
    tthNew = np.arctan(radNew*mmPerPixel/detectDist)
    
    return etaNew, tthNew, deta, dtth

def checkTrack(track_data,scan_ind,scan,frm):
    trackFound = False
    om_ind = 0
    if len(track_data) > scan_ind:
        for om_ind, omTrack in enumerate(track_data[scan_ind]):
            cond1 = omTrack['scan'] == scan
            cond2 = omTrack['frm'] == frm
            if cond1 and cond2:
                trackFound = True
                break
    return trackFound, om_ind

def checkTruth(truth_data,scan_ind,scan,frm):
    truthFound = False
    om_ind = 0
    if len(truth_data[scan_ind]) > scan_ind:
        for om_ind, omTruth in enumerate(truth_data[scan_ind]):
            if len(omTruth) > 0:
                cond1 = omTruth['scan'] == scan
                cond2 = omTruth['frm'] == frm
                if cond1 and cond2:
                    truthFound = True
                    break
    return truthFound, om_ind

def loadROI(dataPath,scan,frm,etaRoi,tthRoi,params):
    if os.path.isdir(dataPath):
        fnames = timeToFile(scan,dataPath)
        isFile = False
    else:
        fnames = dataPath
        isFile = True
        if params['detector'] == 'eiger_sim':
            fnames = fnames.replace('*','{scan}').format(scan=scan)
            roi = loadPolarROI(fnames,tthRoi,etaRoi,frm,params)
            return roi
        
    # Show roi
    if isFile:
        template = dataPath.format(num2=scan)
        fnames = glob.glob(template)
        roi = loadPolarROI(fnames[0],tthRoi,etaRoi,frm,params)
    else:
        roi = loadPolarROI(fnames,tthRoi,etaRoi,frm,params)        
    return roi