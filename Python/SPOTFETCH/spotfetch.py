# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 13:22:47 2024

@author: dpqb1
"""
import numpy as np
import h5py
import hdf5plugin
import warnings
import pickle
import pandas as pd
import os
import time
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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
    tth = kwargs.get('tth', None)
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
    if tth != None:
        tths = spotData['tths']
        cond3 = (tths > tth-dtth) & (tths < tth+dtth)
    else:
        cond3 = True
    if frm != None:
        frms = spotData['ome_idxs']
        cond4 = frms == frm
    else:
        cond4 = True
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
        roiOmega[i] = np.sum(track[i]['roi'],axis=None)
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
        roiOmega[i] = np.sum(track[i]['roi'],axis=None)
        omegaRange[i] = frameToOmega(track[i]['frm'])
        meanOmega += i*roiOmega[i]
    meanOmega = meanOmega/np.sum(roiOmega)
    
    for i in range(len(track)):
        varOmega += roiOmega[i]*(i-meanOmega)**2/np.sum(roiOmega,None)*step**2

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

        # 2. Load ROI and interpolate to polar coordinates 
        if detectFrame == []:
            detectFrame = frame
        if params['detector'] == 'dexela':
            roi_polar = cd.loadDexPolarRoi(fnames,tth,eta,detectFrame,params)
        elif params['detector'] == 'eiger':
            roi_polar = cd.loadEigerPolarRoi(fnames[0],tth,eta,detectFrame,params)
        
        roi_list.append(roi_polar)
        
    return roi_list

def plotROIs(roi_list,num_cols = 5):
    # Create subplots
    num_images = len(roi_list)
    num_rows = (num_images + num_cols - 1) // num_cols  # Calculate number of rows needed

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(10, 10))

    # Plot each image in a subplot
    for i, ax in enumerate(axes.flat):
        if i < num_images:
            ax.imshow(roi_list[i])
            ax.axis('off')  # Turn off axis
            ax.set_title(f'Spot{i+1}')  # Set title for each subplot
            # Add colorbars
            pcm = ax.pcolormesh(roi_list[i])
            fig.colorbar(pcm,ax=ax)
                    
    # Adjust layout to prevent overlapping
    plt.tight_layout()
    
    plt.show()
    
    return fig
    
def plotSpotWedges(spotData,fnames,frame,params,grains=[],detectFrame=[]):
    tths = spotData['tths']
    etas = spotData['etas']     
    ome_idxs = spotData['ome_idxs'] 
    # Get spot indices at frame
    if grains == []:
        spotInds = np.where(ome_idxs == frame)[0]
    else:
        spotInds = findSpots(spotData,grains=grains,frm=frame)
    
    roiSize = params['roiSize']
    imSize = params['imSize']
    center = (imSize[0]//2,imSize[1]//2,)

    if os.path.isfile(fnames):
        fnames = [fnames]
    if detectFrame == []:
        detectFrame = frame
    if params['detector'] == 'dexela':
        b = cd.load_dex(fnames,params,detectFrame)
    elif params['detector'] == 'eiger':
        b = cd.load_eiger(fnames,params,detectFrame)
    
    
    fig = plt.figure()
    # b[b>100] = 100
    plt.imshow(b)
    plt.clim(np.median(b),np.median(b)+10)
    plt.colorbar()
    
    for ind in spotInds:
        # 0. Load spot information
        tth = tths[ind]
        eta = -etas[ind]

        # 1. Load YAML data
        detectDist, mmPerPixel, ff_trans = cd.loadYamlData(params,tth,eta)
        
        # 2. Construct rad, eta domain
        rad_dom, eta_dom = cd.polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
   
        rad = np.round(detectDist*np.tan(tth)/mmPerPixel)
        
        wedge = patches.Wedge([center[1],center[0]],rad+roiSize[0]/2,\
                180/np.pi*eta_dom[0],180/np.pi*eta_dom[-1],\
                linewidth=1,width=roiSize[1],fill=0,color='r')
            
        plt.gca().add_patch(wedge)
        # x = rad*np.cos(eta) + center[1]
        # y = rad*np.sin(eta) + center[0]
        # plt.plot(x,y,color='r',marker='x')
        
    return fig

def roiAdjacent(ind,tth,eta,frame,omFrms,timeFrms,params,dataDir):
    
    fig, axes = plt.subplots(len(omFrms), len(timeFrms), figsize=(len(omFrms), len(timeFrms)))
    
    # add enumerate HERR ND THEN coordinate subplots
    for i,frm in enumerate(omFrms):
        frm = wrapFrame(frm)
        
        for j,t in enumerate(timeFrms):
            fnames = timeToFile(t,dataDir)
            roi = cd.loadPolarROI(fnames,tth,eta,frm,params)
            img = axes[i,j].imshow(roi)
            fig.colorbar(img, ax=axes[i,j])
            axes[i,j].set_title(f'$\omega={frm}$, t={t}')
            
    
def evaluateROI(fnames,prevTracks,tth,eta,frm,scan,params):
    # 0. Parameters
    roiSize = params['roiSize']
    
    # 1. Load ROI
    roi = cd.loadPolarROI(fnames,tth,eta,frm,params)
    tth_vals, eta_vals = np.indices(roi.shape)
    
    # 2. Estimate peak parameters (use from previous timestep)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p0 = fitpeak.estimate_pk_parms_2d(eta_vals,tth_vals,roi,"gaussian")
            p = fitpeak.fit_pk_parms_2d(p0,eta_vals,tth_vals,roi,"gaussian")
            if (p[3]==0): p[3] += 0.001
            if (p[4]==0): p[4] += 0.001
            
    except:
        peakFound = False
        return 0, peakFound
    
    if (p[1] > roiSize[0]-0.5) | (p[2] > roiSize[1]-0.5) | (p[1] < -0.5) | (p[2] < -0.5):
        peakFound = False
        # print('Mean exceeds ROI')
        return 0, peakFound
    
    residual = fitpeak.fit_pk_obj_2d(p,eta_vals,tth_vals,roi,"gaussian")
    rel_error = np.linalg.norm(residual)/np.linalg.norm(roi.flatten())
    
    detectDist, mmPerPixel, ff_trans = cd.loadYamlData(params,tth,eta)
    rad_dom, eta_dom = cd.polarDomain(detectDist,mmPerPixel,tth,eta,roiSize)  
    
    deta = abs(eta_dom[1] - eta_dom[0])
    hypot = detectDist*np.cos(tth)
    dtth = np.arctan(mmPerPixel/hypot)
    
    etaNew = eta_dom[int(np.round(p[1]))]
    radNew = rad_dom[int(np.round(p[2]))]
    tthNew = np.arctan(radNew*mmPerPixel/detectDist)
    
    newTrack = {}
    newTrack['p'] = p
    newTrack['err'] =  rel_error 
    newTrack['etaRoi'] = eta
    newTrack['tthRoi'] = tth  
    newTrack['eta'] =  etaNew
    newTrack['tth'] =  tthNew
    newTrack['frm'] =  frm
    newTrack['scan'] = scan
    newTrack['roi'] = roi
    newTrack['deta'] = deta
    newTrack['dtth'] = dtth
    
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

def scatterOmeg(trackData):
    T = len(trackData)
    K = len(trackData[0])
    
    fig, axes = plt.subplots(K, 1, figsize=(K, 1))
    fig.subplots_adjust(hspace=0.5) 
    
    for t in range(T):
        if trackData[t] != []:
            for k in range(K):
                L = len(trackData[t][k])
                if L > 0:
                    for j in range(L):
                        omega = frameToOmega(trackData[t][k][j]['frm'])
                        axes[k].scatter(t,omega,marker='s',color='b')
                        axes[k].set_title(f'spot {k}')
                        axes[k].set_xlim((0,20))

    return fig

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
    try:
        pool = Pool(params['pool'])
        print(f'{pool._processes} workers')
        parallelFlag = True
    except:
        print('Not parallel')
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
        # if read is successful...
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
            nTrack, peakFound = evaluateROI(fnames,newTrack,\
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
            nTrack, peakFound = evaluateROI(fnames,newTrack,\
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
    
def processSpot(k,t,params,outPath,fnames):
    print(f'{k}', end=" ")
    outFile = os.path.join(outPath,'outputs',f'trackData_{k}.pkl')
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
            try:
                if params['detector'] == 'eiger':
                    ims[frm,:,:]
                    notLoaded = False                    
            except:
                ims = imageseries.open(fnames[0], format='eiger-stream-v1')
                
        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fnames,prevTracks,\
                            tth,eta,int(frm),t,params)
        # Add to list if peakFound
        if peakFound: 
            # print(f'Peak found at frame {frm}')
            trackData[T-1].append(newTrack)
            compareTrack = trackData[T-1]
            frm1 = trackData[T-1][0]['frm']
            frm2 = trackData[T-1][-1]['frm']
            break
    
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
                compareTrack = trackData[T-1]
                frm1 = trackData[T-1][0]['frm']
                frm2 = trackData[T-1][-1]['frm']
                break
    else:
        peakFound = False
    
    # Incremental Search if we have a peak found
    # Search down
    if len(trackData[T-1]) > 0: peakFound = True
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
            
            # Follow up by searching adjacent omegaa for a track
            
    
    try: pool.close()
    except: print('Pool close didnt work')

def compAvgParams(track):   
    J = len(track) 
    avgFWHMeta = 0
    avgFWHMtth = 0
    avgEta = 0
    avgTth = 0
    for j in range(J):
        avgFWHMeta += track[j]['p'][3]/J*track[j]['deta']
        avgFWHMtth += track[j]['p'][4]/J*track[j]['dtth']
        avgEta += track[j]['eta']/J
        avgTth += track[j]['tth']/J
    
    avgEta = avgEta*(180/np.pi)
    avgTth = avgTth*(180/np.pi)
    avgFWHMeta = avgFWHMeta*(180/np.pi)
    avgFWHMtth = avgFWHMtth*(180/np.pi)
    return avgFWHMeta,avgFWHMtth,avgEta,avgTth

def roiTrackVisual(spotInds,spotData,dome,scanRange,trackPath,dataPath,params):
    initData = {
        'tths': spotData['tths'],
        'etas': spotData['etas'],
        'frms': spotData['ome_idxs']
    }
    
    # Choose spot
    for k in spotInds: #range(10):
        print(f'Showing Spot {k}')
        eta0 = initData['etas'][k]
        tth0 = initData['tths'][k]
        frm0 = initData['frms'][k]
        
        # Define Omega range
        om0 = frm0
        omMin = om0-dome
        omMax = om0+dome
        omRange = np.arange(omMin,omMax+1)
        for i,om in enumerate(omRange):
            omRange[i] = wrapFrame(om)
        
        # Define the number of rows and columns
        rows = len(omRange)
        columns = len(scanRange)
        
        # Create a figure and a set of subplots
        fig, axes = plt.subplots(rows, columns, figsize=(20, 15))
        
        # Remove x and y ticks for clarity
        for ax_row in axes:
            for ax in ax_row:
                ax.set_xticks([])
                ax.set_yticks([])
        
        # Add common labels
        fig.text(0.04, 0.5, r'$\omega$ frame', va='center', rotation='vertical', fontsize=24)
        fig.text(0.5, 0.04, 'Scan #', ha='center', fontsize=24)
        fig.text(0.5, 0.95, f'Spot {k}', ha='center', fontsize=32)
        
        # Path to the track data file
        track_file = os.path.join(trackPath,f'trackData_{k}.pkl')
        
        # Initial plot setup
        if os.path.exists(track_file):
            with open(track_file, 'rb') as f:
                track_data = pickle.load(f)
                
        # Organize all tracks
        T = len(track_data)
        omTrack = np.zeros((rows,columns)) - 1
        scanTrack = np.zeros((rows,columns)) - 1
        etaRoiTrack = np.zeros((rows,columns)) + eta0
        tthRoiTrack = np.zeros((rows,columns)) + tth0
        etaTrack = np.zeros((rows,columns)) + eta0
        tthTrack = np.zeros((rows,columns)) + tth0
        p1Track = np.zeros((rows,columns))
        p2Track = np.zeros((rows,columns))
        
        for t in range(T):
            track = track_data[t]
            if len(track) > 0:
                for j in range(len(track)):
                    scan = track[j]['scan']
                    frm = track[j]['frm']
                    
                    if (frm in omRange) & (scan in scanRange):
                        ind1 = int(np.where(scanRange == scan)[0][0])
                        ind2 = int(np.where(omRange == frm)[0][0])
                    
                        eta = track[j]['eta']
                        tth = track[j]['tth']
                        etaRoi = track[j]['etaRoi']
                        tthRoi = track[j]['tthRoi']
                        p = track[j]['p']
                        
                        scanTrack[ind2,ind1] = scan
                        omTrack[ind2,ind1] = frm
                        etaTrack[ind2,ind1] = eta
                        tthTrack[ind2,ind1] = tth
                        etaRoiTrack[ind2,ind1] = etaRoi
                        tthRoiTrack[ind2,ind1] = tthRoi
                        p1Track[ind2,ind1] = p[1]
                        p2Track[ind2,ind1] = p[2]
                    
        for j in range(columns):
            if os.path.isdir(dataPath):
                fnames = timeToFile(scanRange[j],dataPath)
                isFile = False
            else:
                fnames = dataPath
                isFile = True
                
            for i in range(rows): 
                ax = axes[i, j]
                ax.clear()
    
                etaRoi = etaRoiTrack[i,j]
                tthRoi = tthRoiTrack[i,j]
                frm = omRange[i]
                scan = scanRange[j]
                p1 = p1Track[i,j]
                p2 = p2Track[i,j]
                
                # Show roi
                if isFile:
                    print(dataPath)
                    template = dataPath.format(num2=scan)
                    print(template)
                    fnames = glob.glob(template)
                    print(fnames)
                    roi = cd.loadPolarROI(fnames,tthRoi,etaRoi,frm,params)
                else:
                    roi = cd.loadPolarROI(fnames,tthRoi,etaRoi,frm,params)
                ax.imshow(roi)
                
                if (p1 > 0) & (p2 > 0):
                    # Plot rect
                    eta = etaTrack[i,j]
                    tth = tthTrack[i,j]
                    detectDist, mmPerPixel, ff_trans = cd.loadYamlData(params,tth,eta)
                    boxSize = 10
                    start_col = round(p1 - boxSize//2)
                    start_row = round(p2 - boxSize//2)
                    rect = plt.Rectangle((start_col, start_row), boxSize, boxSize,
                                          linewidth=1, edgecolor='r', facecolor='none')
                    ax.add_patch(rect)
                    
                if i == rows-1:
                    ax.set_xlabel(f'{scan}')
                if j == 0:
                    ax.set_ylabel(f'{frm}')
                    
                ax.set_xticks([])
                ax.set_yticks([])
        plt.draw()
        plt.savefig(f'fig_{k}.png')
        plt.close(fig)

        