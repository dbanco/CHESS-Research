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
import h5py
import warnings
import pickle
import pandas as pd
import os
import time
import matplotlib.pyplot as plt
from hexrd.fitting import fitpeak
from hexrd import imageseries
from multiprocessing import Pool
from functools import partial
from pathlib import Path
import glob
import job_manager as jm
import chess_detectors as cd

FRAME1 = 2
NUMFRAMES = 1440
OMEG_RANGE = 360

def omegaToFrame(omega,startFrame=FRAME1,numFrames=NUMFRAMES,omegaRange=OMEG_RANGE,startOmeg = 0):
    step = omegaRange/numFrames
    frame = (np.floor((omega-startOmeg)/step) + startFrame).astype(int)
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
            
        roi_polar = cd.loadPolarROI([fnames],tth,eta,detectFrame,params)    
        
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
     
def assembleTrack(tthRoi,etaRoi,frm,scan,roiSize,p,tth_vals,eta_vals,roi,params):
    
    residual = fitpeak.fit_pk_obj_2d(p,eta_vals,tth_vals,roi,params['peak_func'])
    rel_error = np.linalg.norm(residual)/np.linalg.norm(roi.flatten())

    etaNew, tthNew, deta, dtth = pixToEtaTth(p[1],p[2],tthRoi,etaRoi,params)

    newTrack = {}
    newTrack['p'] = p
    newTrack['err'] =  rel_error 
    newTrack['etaRoi'] = etaRoi
    newTrack['tthRoi'] = tthRoi  
    newTrack['eta'] =  etaNew
    newTrack['tth'] =  tthNew
    newTrack['frm'] =  frm
    newTrack['scan'] = scan
    newTrack['roi'] = roi
    newTrack['deta'] = deta
    newTrack['dtth'] = dtth
    
    return newTrack
    
def evaluateROI(fnames,prevTracks,tthRoi,etaRoi,frm,scan,params):
    # 0. Parameters
    roiSize = params['roiSize']
    
    # 1. Load ROI
    roi = cd.loadPolarROI(fnames,tthRoi,etaRoi,frm,params)
    tth_vals, eta_vals = np.indices(roi.shape)
    
    # 2. Estimate peak parameters (use from previous timestep)
    p, peakFound = fitModel(eta_vals,tth_vals,roi,params)
    
    if peakFound == False:
        return 0, peakFound
   
    # 3. Assemble potential track
    newTrack = assembleTrack(tthRoi,etaRoi,frm,scan,roiSize,p,tth_vals,eta_vals,roi,params)
    
    # 4. Evaluate potential track    
    peakFound = peakDetected(newTrack,prevTracks,params)
    
    # recon = pkfuncs.gaussian2d(p, eta_vals, tth_vals)
    # plt.figure()
    # plt.subplot(2,1,1)
    # plt.imshow(roi)
    # plt.subplot(2,1,2)
    # plt.imshow(recon)
    
    return newTrack, peakFound

def peakDetected(newTrack,prevTracks,params):
    roiSize = params['roiSize']
    p = newTrack['p']
    eta = newTrack['eta']
    tth = newTrack['tth']
    gamma = params['gamma']
    detectDist, mmPerPixel, ff_trans = cd.loadYamlData(params,tth,eta)
    rad_dom, eta_dom = cd.polarDomain(detectDist,mmPerPixel,tth,eta,roiSize)  
    deta = eta_dom[1]-eta_dom[0]
    dtth = hypot = detectDist*np.cos(tth)
    dtth = np.arctan(mmPerPixel/hypot)
    
    if (p[1] > roiSize[0]) | (p[2] > roiSize[1]):
        peakFound = False
        # print('Mean exceeds ROI')
        return peakFound
    
    if len(prevTracks) == 0:
        peakFound = True
        return peakFound
    
    for pTrack in prevTracks:
        
        # Prev track
        pPrev = pTrack['p']
        etaPrev = pTrack['eta']
        tthPrev = pTrack['tth']
        # New track criterion
        crit1 = (abs(eta-etaPrev)/deta < gamma[0]) & (abs(tth-tthPrev)/dtth < gamma[1])
        crit2 = (abs(p[3] - pPrev[3]) < gamma[2]) & (abs(p[4] - pPrev[4]) < gamma[3])
        if crit1 & crit2:
            peakFound = True
            return peakFound
        else:
            peakFound = False
            # print('Track differs from previous')
            # print(f'shiftDist={shiftDist}, sizeDiff={sizeDiff}')
            
    return peakFound

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

def collectSpotsData(outPath,spotsPath):
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
                    id_nums = np.array(df['# ID'])
                    grain_nums = int(file_name[-7:-4])*np.ones(df['pred X'].shape)
                    tths = np.array(df['meas tth'])
                    etas = np.array(df['meas eta'])
                    omes = np.array(df['meas ome'])*180/np.pi
                    created = 1
                else:
                    Xs = np.append(Xs, np.array(df['pred X']))
                    Ys = np.append(Ys, np.array(df['pred Y']))
                    id_nums = np.append(id_nums, np.array(df['# ID']))
                    grain_nums = np.append(grain_nums,int(file_name[-7:-4])*np.ones(df['pred X'].shape))
                    tths = np.append(tths, np.array(df['meas tth']))
                    etas = np.append(etas, np.array(df['meas eta']))
                    omes = np.append(omes, np.array(df['meas ome'])*180/np.pi)
                
    invalid_nums = (id_nums == -999) | np.isnan(tths.astype(float)) | np.isnan(etas.astype(float))       
    Xs = np.delete(Xs, invalid_nums).astype(float)
    Ys = np.delete(Ys, invalid_nums).astype(float)
    id_nums = np.delete(id_nums, invalid_nums).astype(float)
    grain_nums = np.delete(grain_nums, invalid_nums).astype(float)
    tths = np.delete(tths, invalid_nums).astype(float)
    etas = np.delete(etas, invalid_nums).astype(float)
    omes = np.delete(omes, invalid_nums).astype(float)
    
    ome_idxs = omegaToFrame(omes)
    
    sort_ind = np.argsort(omes)
    
    ome_idxs = ome_idxs[sort_ind]
    Xs = Xs[sort_ind]
    Ys = Ys[sort_ind]
    id_nums = id_nums[sort_ind]
    tths = tths[sort_ind]
    etas = etas[sort_ind]
    omes = omes[sort_ind]
    grain_nums = grain_nums[sort_ind]
    saveFile = os.path.join(outPath,'spots.npz')
    np.savez(saveFile,Xs=Xs,Ys=Ys,id_nums=id_nums,\
    tths=tths,etas=etas,omes=omes,ome_idxs=ome_idxs,grain_nums=grain_nums)

def spotTracker(dataPath,outPath,spotData,spotInds,params,num1,num2,advance):
    # i = 1
    # t = scan1
    newData = False
    if params['parallelFlag']:
        try:
            pool = Pool(params['pool'])
            print(f'{pool._processes} workers')
            parallelFlag = True
        except:
            print('Not parallel')
            parallelFlag = False
    else:
        parallelFlag = False
    
    while True:
        # Try reading in file for new scan  
        try:
            # if os.path.isfile(dataPath):
            template = dataPath.format(num2=num2)
            template2 = dataPath.format(num2=num2+1)
            pattern = os.path.join(template)
            pattern2 = os.path.join(template2)

            # Use glob to find files that match the pattern
            fnames = glob.glob(pattern)
            fnames2 = glob.glob(pattern2)
            
            if fnames2 == []:
                time.sleep(1)
                continue
            elif Path(fnames2[0]).stat().st_size >= 2e6:
                print('File ready')
            else:
                time.sleep(1)
                continue
            # else:
                # dataDir = os.path.join(dataPath,f'{num2}','ff')
                # fnames = pathToFile(dataDir)
        except:
            if newData:
                print('No new data',end=" ")
                newData = False
            else:
                print('.',end="")
            time.sleep(1)
            continue
        # if the files exist
        if len(fnames) > 0:
            newData = True        
            
            print(f'Scan {num2}, Spot:', end=" ")
            if parallelFlag:
                pool.starmap(partial(processSpot,t=num2,params=params,\
                                      outPath=outPath,fnames=fnames),zip(spotInds))
            else:
                for s,k in enumerate(spotInds): 
                    processSpot(k,num2,params,outPath,fnames)
            if advance:
                num1 += 1
                num2 += 1
            else:
                break
        
def initSpot(k,etaRoi,tthRoi,frm,t,params,outPath,fnames):
    print(k)
    prevTracks = []
    frm = int(frm)
    newTrack, peakFound = evaluateROI(fnames,prevTracks,tthRoi,etaRoi,int(frm),t,params)
    trackData = []
    if peakFound:
        trackData.append([newTrack])
        
        if True:#numFrames == NUMFRAMES:
             wrap = True
        else:
             wrap = False
        
        # Search down
        frm1 = frm
        while peakFound:
            frm1 = frm1 - 1
            if wrap:
                frm = int(wrapFrame(frm1))
            else:
                frm = frm1
                if frm < 0: break
            # Load ROI and fit peak
            nTrack, peakFound = evaluateROI(fnames,[newTrack],\
                                              tthRoi,etaRoi,frm,t,params)
            # Add to list if peakFound
            if peakFound: 
                # print(f'Found more at {frm1}')
                trackData[0].insert(0,nTrack)
                if len(trackData[0]) > 10:
                    break
    
        # Search up
        frm2 = frm
        while peakFound:
            frm2 = frm2 + 1
            frm = int(wrapFrame(frm2))
            if wrap:
                frm = int(wrapFrame(frm2))
            else:
                frm = frm2
                if frm > NUMFRAMES-1: break
            # Load ROI and fit peak
            nTrack, peakFound = evaluateROI(fnames,[newTrack],\
                                              tthRoi,etaRoi,frm,t,params)
            # Add to list if peakFound
            if peakFound: 
                # print(f'Found more at {frm2}')
                trackData[0].append(nTrack)
                if len(trackData[0]) > 10:
                    break
        
    outFile = os.path.join(outPath,f'trackData_{k}.pkl')
    with open(outFile, 'wb') as f:
        pickle.dump(trackData, f)
    
def processSpot(k,t,params,trackPath,fnames):
    print(f'{k}', end=" ")
    outFile = os.path.join(trackPath,f'trackData_{k}.pkl')
    # Load in track so far
    with open(outFile, 'rb') as f:
        trackData = pickle.load(f)
    
    if len(trackData) == 1 and trackData[0] == []:
        return 0
    
    trackData.append([])
    T = len(trackData)
    
    # Determine if scan wraps in omega or if LODI
    if params['detector'] == 'eiger':
        ims = imageseries.open(fnames[0], format='eiger-stream-v1')
        numFrames = ims.shape[0]
    elif params['detector'] == 'dexela':
        with h5py.File(fnames[0], 'r') as file1:
            numFrames = file1['/imageseries/images'].shape[0]
    elif params['detector'] == 'eiger_sim':
        simData = np.load(fnames[0])
        numFrames = simData['nframes']
                        
    if True:#numFrames == NUMFRAMES:
        wrap = True
    else:
        wrap = False

    prevTracks = []
    lag = 1
    while (len(prevTracks) == 0) & (T-lag >= 0):
        prevTracks = trackData[T-lag]
        lag = lag + 1
        
    # Initial Search: through all current omega tracks, then check up and down for\
    # tracks (not sure exactly when search through omega will be considered done)    
    for track in prevTracks:
        eta = track['eta']
        tth = track['tth']
        frm = track['frm']
        
        notLoaded = True
        while notLoaded:
            if params['detector'] == 'eiger':
                try:
                    ims[frm,:,:]
                    notLoaded = False                    
                except:
                        ims = imageseries.open(fnames[0], format='eiger-stream-v1')
            elif params['detector'] == 'eiger_sim':
                try:
                    simData[f'{frm}_row']
                    notLoaded = False                    
                except:
                        simData = np.load(fnames[0])
                
        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fnames,prevTracks,\
                            tth,eta,int(frm),t,params)
        # Add to list if peakFound
        if peakFound: 
            # print(f'Peak found at frame {frm}')
            trackData[T-1].append(newTrack)
            frm1 = trackData[T-1][0]['frm']
            frm2 = trackData[T-1][-1]['frm']
            
    
    # Conduct Expanded Search if no peaks were found
    if len(trackData[T-1]) == 0 & len(prevTracks) > 0:
        frm1 = prevTracks[0]['frm']
        frm2 = prevTracks[-1]['frm']
        expandRange = list(range(frm1-3,frm1)) + list(range(frm2+1,frm2+4))
        for frm in expandRange:
            if wrap:
                frm = int(wrapFrame(frm))
            else:
                if frm < 0 | frm > NUMFRAMES-1: break
            newTrack, peakFound = evaluateROI(fnames,prevTracks,\
                                tth,eta,frm,t,params)
            if peakFound: 
                # print(f'Peak found at frame {frm}')
                trackData[T-1].append(newTrack)
                frm1 = trackData[T-1][0]['frm']
                frm2 = trackData[T-1][-1]['frm']
                break
    else:
        peakFound = False
    
    # Incremental Search if we have a peak found
    # Search down
    if len(trackData[T-1]) > 0: 
        peakFound = True
        compareTrack = trackData[T-1]
    while peakFound:
        frm1 = frm1 - 1
        if wrap:
            frm = int(wrapFrame(frm1))
        else:
            frm = frm1
            if frm < 0: break
        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fnames,compareTrack,\
                                          tth,eta,frm,t,params)
        # Add to list if peakFound
        if peakFound: 
            # print(f'Found more at {frm1}')
            trackData[T-1].insert(0,newTrack)
            if len(trackData[T-1]) > 10:
                break

    # Search up
    if len(trackData[T-1]) > 0: peakFound = True
    while peakFound:
        frm2 = frm2 + 1
        frm = int(wrapFrame(frm2))
        if wrap:
            frm = int(wrapFrame(frm2))
        else:
            frm = frm2
            if frm > NUMFRAMES-1: break
        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fnames,compareTrack,\
                                          tth,eta,frm,t,params)
        # Add to list if peakFound
        if peakFound: 
            # print(f'Found more at {frm2}')
            trackData[T-1].append(newTrack)
            if len(trackData[T-1]) > 10:
                break

    with open(outFile, 'wb') as f:
        pickle.dump(trackData, f)
        
def processSpot_wrapper(inputFile):
    with open(inputFile, 'rb') as f:
        k,t,params,outPath,fnames = pickle.load(f) 
    processSpot(k,t,params,outPath,fnames)    
        
def spotTrackerJobs(dataPath, outPath, exsituPath, spotData, spotInds, params, scan1):
    # Job template
    job_script_template = """#!/bin/bash
#$ -N spotTrack_{k}
#$ -cwd
#$ -l h_vmem=4G
#$ -l h_rt=1:00:00
#$ -j y
#$ -o spotTrack_{k}.out
source hexrdenv/bin/activate
conda init bash
conda activate hexrd-env
python3 -c "import sys; sys.path.append('CHESS-Research/Python/SPOTFETCH/'); import spotfetch as sf; sf.processSpotJob('{inputFile}')"
"""

    # Initialize 
    initData = {}
    initData['tths'] = spotData['tths'][spotInds]
    initData['etas'] = spotData['etas'][spotInds]
    initData['frms'] = spotData['ome_idxs'][spotInds]
    
    # Initialize tracks
    fnames = pathToFile(exsituPath)

    i = 0
    t = scan1-1  # Scan index
    print('')
    print(f'Scan {t}, Spot:', end=" ")
    for s, k in enumerate(spotInds):
        print(f'{k}', end=" ")
        etaRoi = initData['etas'][s]
        tthRoi = initData['tths'][s]
        frm = initData['frms'][s]
        initSpot(k, etaRoi, tthRoi, frm,t, params, outPath, fnames)
    
    i += 1
    t += 1
    while True:
        # Try reading in file for new scan
        dataDir = os.path.join(dataPath,f'{t}','ff')
        try:
            fnames = pathToFile(dataDir)
        except:
            print('No new data')
            time.sleep(1)
            continue
        # If read is successful...
        if len(fnames) > 0:

            print('')
            print(f'Scan {t}, Spot:', end=" ")
            for s, k in enumerate(spotInds): 
                print(f'{k}', end=" ")
                # Save input file
                inputFile = os.path.join(outPath,'inputs',f'inputs_{k}.pkl')
                with open(inputFile, 'wb') as f:
                    pickle.dump((k,s,t,params,outPath,fnames), f)
                # Write job .sh file
                job_script = job_script_template.format(k=k,inputFile=inputFile)
                script_filename = os.path.join(outPath,'jobs',f'job_{k}.sh')
                with open(script_filename, 'w') as f:
                    f.write(job_script)
                # Submit job
                os.system(f"qsub {script_filename}")
 
            # Wait for my jobs to finish
            job_ids = jm.collect_job_ids()
            jm.wait_for_jobs(job_ids)
                
            i += 1
            t += 1
            

def initExsituTracks(outPath,exsituPath,spotData,spotInds,params,scan0):
    # Initialize 
    initData = {}
    initData['tths'] = spotData['tths'][spotInds]
    initData['etas'] = spotData['etas'][spotInds]
    initData['frms'] = spotData['ome_idxs'][spotInds]
    
    # Initialize tracks
    if os.path.isfile(exsituPath):
        fnames = [exsituPath]
    else:
        fnames = pathToFile(exsituPath)
    
    if params['parallelFlag']:
        try:
            pool = Pool(params['pool'])
            print(f'{pool._processes} workers, Ex-Situ Scan ({scan0})')
            inVars = zip(spotInds,initData['etas'],initData['tths'],initData['frms'])
            pool.starmap(partial(initSpot,t=scan0,params=params,outPath=outPath,\
                          fnames=fnames),inVars)
        except:
            # Scan index
            print('Not parallelized')
            print(f'Ex-Situ Scan ({scan0}), Spot:', end=" ")
            for s,k in enumerate(spotInds):
                print(f'{k}', end=" ")
                etaRoi = initData['etas'][s]
                tthRoi = initData['tths'][s]
                frm = initData['frms'][s]
                initSpot(k,etaRoi,tthRoi,frm,scan0,params,outPath,fnames)
                
        try: pool.close()
        except: print('Pool close didnt work')
    else:
        # Scan index
        print('Not parallelized')
        print(f'Ex-Situ Scan ({scan0}), Spot:', end=" ")
        for s,k in enumerate(spotInds):
            print(f'{k}', end=" ")
            etaRoi = initData['etas'][s]
            tthRoi = initData['tths'][s]
            frm = initData['frms'][s]
            initSpot(k,etaRoi,tthRoi,frm,scan0,params,outPath,fnames)
            
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


        
def trackingResults(spotInds,spotFiles,scanRange,trackPath,dataPath,params):
    state0Data = np.load(spotFiles[0])
    initData = {
        'tths': state0Data['tths'],
        'etas': state0Data['etas'],
        'frms': state0Data['ome_idxs']
    }
    
    numSpots =  len(spotInds)
    T = len(scanRange)
    foundInd = np.zeros((numSpots,T))
    truthDist = np.zeros((numSpots,T))
    changeDist = np.zeros((numSpots,T))
    
    foundInd[:] = np.nan
    truthDist[:] = np.nan
    changeDist[:] = np.nan
    
    spotData = []
    for i in range(len(spotFiles)):
        spotData.append(np.load(spotFiles[i]))
        
    for k_ind,k in enumerate(spotInds): #range(10):
        print(f'Spot {k}')
        eta0 = initData['etas'][k]
        tth0 = initData['tths'][k]
        frm0 = initData['frms'][k]
        
        # Load track data for spot
        track_file = os.path.join(trackPath,f'trackData_{k}.pkl')
        if os.path.exists(track_file):
            with open(track_file, 'rb') as f:
                track_data = pickle.load(f)
        
        # Metrics to compute: (Need to gather spot data)
        #   Tracks found at steps 0-4 (binary array for each spot)
        #   Compute distance of found spot from spot.out files
        #   Compute distance from init spot location of spot.out file
        
        # 1. Identify # of omega tracks per time step
        if len(track_data) > 0:
            for track in track_data:
                if len(track) > 0:
                    # 1. Count omegas at each spot/time
                    foundInd[k_ind,track[0]['scan']] = len(track)
                    # 2. Compute change in eta/tth over time
                    etaChange = track[0]['eta']-track[0]['etaRoi']
                    tthChange = track[0]['tth']-track[0]['tthRoi']
                    changeDist[k_ind,track[0]['scan']] = np.sqrt(etaChange**2 + tthChange**2)                    
                    # 3. Compute difference from true positions
                    etaTruth = track[0]['eta'] - spotData[track[0]['scan']]['etas'][k_ind]
                    tthTruth = track[0]['tth'] - spotData[track[0]['scan']]['tths'][k_ind]                  
                    truthDist[k_ind,track[0]['scan']] = np.sqrt(etaTruth**2 + tthTruth**2)
                    # 4. Count number of spots with tracks at initial time step                 
                    numTracks0 = np.sum(foundInd[:,0] > 0)
    return foundInd, numTracks0, changeDist, truthDist         
       
def loadROI(dataPath,scan,frm,etaRoi,tthRoi,params):
    if os.path.isdir(dataPath):
        fnames = timeToFile(scan,dataPath)
        isFile = False
    else:
        fnames = dataPath
        isFile = True
        if params['detector'] == 'eiger_sim':
            fnames = fnames.replace('*','{scan}').format(scan=scan)
            roi = cd.loadPolarROI(fnames,tthRoi,etaRoi,frm,params)
            return roi
        
    # Show roi
    if isFile:
        template = dataPath.format(num2=scan)
        fnames = glob.glob(template)
        roi = cd.loadPolarROI(fnames,tthRoi,etaRoi,frm,params)
    else:
        roi = cd.loadPolarROI(fnames,tthRoi,etaRoi,frm,params)        
    return roi

def showROI(ax,dataPath,scan,frm,tthRoi,etaRoi,params):
    
    roi = loadROI(dataPath,scan,frm,etaRoi,tthRoi,params)
        
    detectDist, mmPerPixel, ff_trans, ff_tilt = cd.loadYamlData(params,tthRoi,etaRoi)
    rad_dom, eta_dom = cd.polarDomain(detectDist, mmPerPixel,\
                          tthRoi, etaRoi, params['roiSize'])
    tth_dom = np.arctan(rad_dom*mmPerPixel/detectDist)
    
    ax.imshow(roi)
    plt.xticks([0,39], [f'{eta_dom[0]:.4g}',f'{eta_dom[-1]:.4g}'])
    plt.yticks([0,39], [f'{tth_dom[0]:.4g}',f'{tth_dom[-1]:.4g}'])
  
    ax.text(1, 2, f'max: {roi.max():.2f}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 4, f'scan: {scan}', color='white', fontsize=12, weight='bold')
    # ax.text(1, 6, f'frame: {frm}', color='white', fontsize=12, weight='bold')
    ax.text(1, 4, f'eta: {etaRoi:.4f}', color='white', fontsize=12, weight='bold')
    ax.text(1, 6, f'tth: {tthRoi:.4f}', color='white', fontsize=12, weight='bold')

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
        ax.plot(x_pos,y_pos,marker='s',markersize=8,\
                fillstyle='none',color='red')
    ax.plot(x_pos,y_pos,marker='x',markersize=8,color='red')
    # Previous track
    [y_pos, x_pos] = etaTthToPix(etaPrev,tthPrev,etaRoi,tthRoi,params)
    if (y_pos >= 0) and (y_pos < params['roiSize'][0]) and\
       (x_pos >= 0) and (x_pos < params['roiSize'][1]):
        ax.plot(x_pos,y_pos,marker='s',markersize=10,\
                fillstyle='none',color='yellow')


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
    
def etaTthToPix(eta,tth,etaRoi,tthRoi,params):
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans, ff_tilt = cd.loadYamlData(params,tthRoi,etaRoi)
    rad_dom, eta_dom = cd.polarDomain(detectDist, mmPerPixel,\
                          tthRoi, etaRoi, roiSize)
        
    rad = np.tan(tth)*detectDist/mmPerPixel
    
    row_pos = (rad-rad_dom[0])/(rad_dom[-1]-rad_dom[0])*(roiSize[0]-1)
    col_pos = (eta-eta_dom[0])/(eta_dom[-1]-eta_dom[0])*(roiSize[1]-1)
    
    return row_pos, col_pos

def pixToEtaTth(p1,p2,tthRoi,etaRoi,params):
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans = cd.loadYamlData(params,tthRoi,etaRoi)
    rad_dom, eta_dom = cd.polarDomain(detectDist,mmPerPixel,tthRoi,etaRoi,roiSize)  
    
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
    if len(track_data) > 0:
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
    if len(truth_data[scan_ind]) > 0:
        for om_ind, omTruth in enumerate(truth_data[scan_ind]):
            if len(omTruth) > 0:
                cond1 = omTruth['scan'] == scan
                cond2 = omTruth['frm'] == frm
                if cond1 and cond2:
                    truthFound = True
                    break
    return truthFound, om_ind

