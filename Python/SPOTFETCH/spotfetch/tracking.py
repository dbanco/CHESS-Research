# -*- coding: utf-8 -*-
"""
spotfetch.tracking

Tools for visualizing X-ray diffraction data and results.

Created on: Fri Nov 15 23:00:28 2024
Author: Daniel Banco
"""
import numpy as np
import h5py
import pickle
import os
import time
from hexrd.fitting import fitpeak
from hexrd import imageseries
from multiprocessing import Pool
from functools import partial
from pathlib import Path
import glob
import job_manager as jm
import spotfetch as sf
from scipy.ndimage import gaussian_filter
from skimage.feature import hessian_matrix
from scipy.ndimage import label

FRAME1 = 2
NUMFRAMES = 1440
OMEG_RANGE = 360

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
            
def stepOmegaFrame(frm1,step,wrap):
    frm1 = frm1 + step
    if wrap:
        frm = int(sf.wrapFrame(frm1))
    else:
        frm = frm1
    return frm

def wrapOmegaCheck(fnames,params):
    
    if params['detector'] == 'dexela':
        with h5py.File(fnames[0], 'r') as file1:
            numFrames = file1['/imageseries/images'].shape[0]
    elif params['detector'] == 'eiger':
        ims = imageseries.open(fnames, format='eiger-stream-v1')
        numFrames = ims.shape[0]
    elif params['detector'] == 'eiger_sim':
        simData = np.load(fnames)
        numFrames = simData['nframes']
                        
    if numFrames == NUMFRAMES:
        wrap = True
    else:
        wrap = False
        
    return wrap
   
def searchDown(peakFound,frm1,wrap,fnames,compareTrack,tth,eta,scan,params,trackData,T):
    while peakFound:
        # Increment Omega frame
        frm1 = stepOmegaFrame(frm1,-1,wrap)
        if frm1 < FRAME1: break

        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fnames,compareTrack,tth,eta,frm1,scan,params,omegaCompare=True)
        # Add to list if peakFound
        if peakFound: 
            # print(f'Found more at {frm1}')
            trackData[T-1].insert(0,newTrack)
            if len(trackData[T-1]) > 10:
                break
    return trackData

def searchUp(peakFound,frm2,wrap,fnames,compareTrack,tth,eta,scan,params,trackData,T,omegaCompare=True):
    while peakFound:
        # Increment Omega frame
        frm2 = stepOmegaFrame(frm2,1,wrap)
        if frm2 > FRAME1+NUMFRAMES-1: break
    
        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fnames,compareTrack,tth,eta,frm2,scan,params,omegaCompare)
        
        # Add to list if peakFound
        if peakFound: 
            # print(f'Found more at {frm2}')
            trackData[T-1].append(newTrack)
            if len(trackData[T-1]) > 10:
                break
    return trackData

def initSpot(spotInd,eta,tth,frm,fnames,trackFile,params):
    if params['benchmarkFlag']:
        start_time = time.time()
    
    outData = {'trackData': [],
               'trackTimes': []}
    prevTracks = [] 
    
    newTrack, peakFound = evaluateROI(fnames,prevTracks,tth,eta,frm,0,params)

    if peakFound:
        outData['trackData'].append([newTrack])
        
        # Determine if scan wraps in omega or if LODI
        wrap = wrapOmegaCheck(fnames,params)
        
        # Search down
        outData['trackData'] = searchDown(peakFound,frm,wrap,fnames,[newTrack],
                                          tth,eta,0,params,outData['trackData'],1)
        # Search up
        outData['trackData'] = searchUp(peakFound,frm,wrap,fnames,[newTrack],
                                        tth,eta,0,params,outData['trackData'],1)
        if params['benchmarkFlag']:
            end_time = time.time()
            outData['trackTimes'].append(end_time-start_time)

    with open(trackFile, 'wb') as f:
        pickle.dump(outData,f)
        
def processSpot(spotInd,t_ind,fnames,trackFile,params):
    # print(f'{spotInd}', end=" ")

    # Load in track so far
    if params['benchmarkFlag']:
        start_time = time.time()
        
    with open(trackFile, 'rb') as f:
        outData = pickle.load(f)
   
    if len(outData['trackData']) == 1 and outData['trackData'][0] == []:
        return 0
    
    outData['trackData'].append([])
    T = len(outData['trackData'])
    
    # Determine if scan wraps in omega or if LODI
    wrap = wrapOmegaCheck(fnames,params)

    prevTracks = []
    lag = 1
    while (len(prevTracks) == 0) & (T-lag >= 0):
        prevTracks = outData['trackData'][T-lag]
        lag = lag + 1
        
    # Initial Search: through all current omega tracks, then check up and down for\
    # tracks (not sure exactly when search through omega will be considered done)    
    peakFound = False
    for track in prevTracks:
        eta = track['eta']
        tth = track['tth']
        frm = track['frm']
        
        # Check that file can be loaded??
        # testLoadFile()        
        
        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fnames,prevTracks,\
                            tth,eta,int(frm),t_ind,params)
            
        # Add to list if peakFound
        if peakFound: 
            # print(f'Peak found at frame {frm}')
            outData['trackData'][T-1].append(newTrack)
            frm1 = outData['trackData'][T-1][0]['frm']
            frm2 = outData['trackData'][T-1][-1]['frm']       
    
    # Conduct Expanded Search if no peaks were found
    if (not peakFound) & len(prevTracks) > 0:
        frm1 = prevTracks[0]['frm']
        frm2 = prevTracks[-1]['frm']
        expandRange = list(range(frm1-3,frm1)) + list(range(frm2+1,frm2+4))
        for frm_i in expandRange:
            if wrap:
                frm_i = int(sf.wrapFrame(frm_i))
            else:
                if frm < 0 | frm > NUMFRAMES-1: break
            newTrack, peakFound = evaluateROI(fnames,prevTracks,\
                                tth,eta,frm_i,t_ind,params)
            if peakFound: 
                # print(f'Peak found at frame {frm}')
                outData['trackData'][T-1].append(newTrack)
                frm1 = outData['trackData'][T-1][0]['frm']
                frm2 = outData['trackData'][T-1][-1]['frm']
                break
    
    # Incremental Search if we have a peak found
    # Search down and up
    if peakFound: 
        compareTrack = outData['trackData'][T-1]
        outData['trackData'] = searchDown(peakFound,frm1,wrap,fnames,compareTrack,
                                          tth,eta,t_ind,params,outData['trackData'],T)
        outData['trackData'] = searchUp(peakFound,frm2,wrap,fnames,compareTrack,
                                        tth,eta,t_ind,params,outData['trackData'],T)

    if params['benchmarkFlag']:
        end_time = time.time()
        outData['trackTimes'].append(end_time-start_time)
        with open(trackFile, 'wb') as f:
            pickle.dump(outData,f)
    else:
        with open(trackFile, 'wb') as f:
            pickle.dump(outData['trackData'],f)
        
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
    fnames = sf.pathToFile(exsituPath)

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
            fnames = sf.pathToFile(dataDir)
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
        fnames = sf.pathToFile(exsituPath)
    
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

def assembleTrack(tthRoi,etaRoi,frm,scan,roiSize,p,residual,roi,params):
    
    rel_error = np.linalg.norm(residual)/np.linalg.norm(roi.flatten())
    etaNew, tthNew, deta, dtth = sf.pixToEtaTth(p[1],p[2],tthRoi,etaRoi,params)

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
    
def evaluateROI(fnames,prevTracks,tth,eta,frm,scan,params,omegaCompare=False):
    # 0. Parameters
    roiSize = params['roiSize']
    
    # 1. Load ROI
    roi = sf.loadPolarROI(fnames,tth,eta,frm,params)
    
    # 2. Estimate peak parameters (use from previous timestep)
    p, peakFound, residual = sf.fitModel(roi,params,tth,eta)
    
    if peakFound == False:
        return 0, peakFound
   
    # 3. Assemble potential track
    newTrack = assembleTrack(tth,eta,frm,scan,roiSize,p,residual,roi,params)
    
    # 4. Evaluate potential track    
    peakFound = peakDetected(newTrack,prevTracks,params,omegaCompare)
    
    return newTrack, peakFound

def peakDetected(newTrack,prevTracks,params,omegaCompare=False):
    p = newTrack['p']
    eta = newTrack['eta']
    tth = newTrack['tth']
    
    if len(prevTracks) == 0:
        peakFound = True
        return peakFound
    
    for pTrack in prevTracks:
        
        criterion = matchCriterion(eta,tth,p,pTrack,params,omegaCompare)
        
        if criterion:
            peakFound = True
            return peakFound
        else:
            peakFound = False
            
    return peakFound

def matchCriterion(eta,tth,p,pTrack,params,omegaCompare=False):
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans, ff_tilt = sf.loadYamlData(params,tth,eta)
    rad_dom, eta_dom = sf.polarDomain(detectDist,mmPerPixel,tth,eta,roiSize)  
    deta = eta_dom[1]-eta_dom[0]
    dtth = hypot = detectDist*np.cos(tth)
    dtth = np.arctan(mmPerPixel/hypot)
    # Prev track
    pPrev = pTrack['p']
    etaPrev = pTrack['eta']
    tthPrev = pTrack['tth']
    
    # New track criterion
    gamma = params['gamma']
    if omegaCompare:
        etaThresh = gamma[4]
        tthThresh = gamma[5]
    else:
        etaThresh = gamma[0]
        tthThresh = gamma[1]
    if params['peak_func'] == 'gaussian':
        crit1 = (abs(eta-etaPrev)/deta < etaThresh) & (abs(tth-tthPrev)/dtth < tthThresh)
        crit2 = (abs(p[3] - pPrev[3]) < gamma[2]) & (abs(p[4] - pPrev[4]) < gamma[3])
        return crit1 & crit2
    elif params['peak_func'] == 'gaussian_rot':
        crit1 = (abs(eta-etaPrev)/deta < etaThresh) & (abs(tth-tthPrev)/dtth < tthThresh)
        maxP = max(p[3],p[4])
        minP = min(p[3],p[4])
        maxPprev = max(pPrev[3],pPrev[4])
        minPprev = min(pPrev[3],pPrev[4])
        crit2 = (abs(minP - minPprev) < gamma[2]) & (abs(maxP - maxPprev) < gamma[3])
        return crit1 & crit2
    

def trackingResults(spotInds,spotFiles,scanRange,trackPath,dataPath,params):   
    numSpots =  len(spotInds)
    T = len(scanRange)
    spotAtTrueFrm = np.zeros((numSpots,T))
    omegaCount = np.zeros((numSpots,T))
    truthDist = np.zeros((numSpots,T))
    changeDist = np.zeros((numSpots,T))
    
    spotAtTrueFrm[:] = False
    omegaCount[:] = np.nan
    truthDist[:] = np.nan
    changeDist[:] = np.nan
    
    spotData = []
    for i in range(len(spotFiles)):
        spotData.append(np.load(spotFiles[i]))
        
    for k_ind,k in enumerate(spotInds): #range(10):
        print(f'Spot {k}')
        x = spotData['Xm'][k]
        y = spotData['Ym'][k]
        eta0, tth0 = sf.xyToEtaTthRecenter(x,y,params)
        
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
                    omegaCount[k_ind,track[0]['scan']] = len(track)
                    # 2. Compute change in eta/tth over time
                    etaChange = track[0]['eta']-track[0]['etaRoi']
                    tthChange = track[0]['tth']-track[0]['tthRoi']
                    changeDist[k_ind,track[0]['scan']] = np.sqrt(etaChange**2 + tthChange**2)                    
                    # 3. Compute difference from true positions
                    etaTruth = track[0]['eta'] - spotData[track[0]['scan']]['etas'][k_ind]
                    tthTruth = track[0]['tth'] - spotData[track[0]['scan']]['tths'][k_ind]                  
                    truthDist[k_ind,track[0]['scan']] = np.sqrt(etaTruth**2 + tthTruth**2)
                    # 4. Count number of spots with tracks at initial time step                 
                    numTracks0 = np.sum(omegaCount[:,0] > 0)
    return omegaCount, numTracks0, changeDist, truthDist     

def trackingResultsSim(spotInds,spotFiles,scanRange,trackPath,params):  
    numSpots =  len(spotInds)
    T = len(scanRange)
    
    spotAtTrueFrm = np.zeros((numSpots,T))
    omegaCount = np.zeros((numSpots,T))
    truthDist = np.zeros((numSpots,T))
    changeDist = np.zeros((numSpots,T))
    trackTimes = np.zeros((numSpots,T))
    
    spotAtTrueFrm[:] = np.nan
    omegaCount[:] = np.nan
    truthDist[:] = np.nan
    changeDist[:] = np.nan
    trackTimes[:] = np.nan
    
    spotData = []
    for i in range(len(spotFiles)):
        spotData.append(np.load(spotFiles[i]))
        
    for k_ind,k in enumerate(spotInds): #range(10):
        print(f'Spot {k}')
        x = spotData[0]['Xm'][k]
        y = spotData[0]['Ym'][k]
        eta0, tth0 = sf.xyToEtaTthRecenter(x,y,params)
        spot_id = sf.getSpotID(spotData[0],k)
                
        # Skip if truth is lost at some point            
        spot_id = sf.getSpotID(spotData[0],k)
        for t in range(1,5):        
            m_ind = sf.matchSpotID(spotData[t],spot_id,k)
            if np.isnan(m_ind):
                skipSpot = True
                break
            else:
                skipSpot = False     
        if skipSpot:
            continue
        
        # Load track data for spot
        track_file = os.path.join(trackPath,f'trackData_{k}.pkl')
        if os.path.exists(track_file):
            with open(track_file, 'rb') as f:
                outData = pickle.load(f)
                track_data = outData['trackData']
                track_times = outData['trackTimes']
        
        # Metrics to compute: (Need to gather spot data)
        #   Tracks found at steps 0-4 (binary array for each spot)
        #   Compute distance of found spot from spot.out files
        #   Compute distance from init spot location of spot.out file      
        
        # Loop over time steps
        if len(track_data) > 0:
            for track in track_data:
                if len(track) > 0:
                    # 0. Get truth information
                    scan_t = track[0]['scan']
                    # 1. Count omegas at each spot/time
                    omegaCount[k_ind,scan_t] = len(track)
                    # 1.2 Get compute times
                    trackTimes[k_ind,scan_t] = track_times[scan_t]
                    # 1.5 Compute the following if spot is in same frame as truth
                    m_ind = sf.matchSpotID(spotData[scan_t],spot_id,k)
                    if not np.isnan(m_ind):
                        x = spotData[scan_t]['Xm'][m_ind]
                        y = spotData[scan_t]['Ym'][m_ind]
                        etaTrue, tthTrue = sf.xyToEtaTthRecenter(x,y,params)
                        frmTrue = spotData[scan_t]['ome_idxs'][m_ind]
                        # 1 Identify which track is at true frame
                        for ind, omTr in enumerate(track):
                            if omTr['frm'] == frmTrue:
                                spotAtTrueFrm[k_ind,scan_t] = True
                                break
                        if np.isnan(spotAtTrueFrm[k_ind,scan_t]):
                            continue
                        # 2. Compute change in eta/tth over time
                        etaChange = track[ind]['eta']-track[ind]['etaRoi']
                        tthChange = track[ind]['tth']-track[ind]['tthRoi']
                        changeDist[k_ind,scan_t] = np.sqrt(etaChange**2 + tthChange**2)  
                        # 3. Compute difference from true positions
                        etaTruth = track[ind]['eta'] - etaTrue
                        tthTruth = track[ind]['tth'] - tthTrue            
                        truthDist[k_ind,scan_t] = np.sqrt(etaTruth**2 + tthTruth**2)
                        
    resultsData = {
    "spotAtTrueFrm": spotAtTrueFrm,
    "omegaCount": omegaCount,
    "trackTimes": trackTimes,
    "changeDist": changeDist,
    "truthDist": truthDist
}
    
    return resultsData

def trackSpot(spotInd,spotData,dataFileSequence,trackPath,params):
    
    trackFile = os.path.join(trackPath,f'trackData_{spotInd}.pkl')
    if os.path.exists(trackFile):
        print('File exists. Spot {spotInd} not tracked.')
        return 0
    
    # Initial spot location 
    x = spotData['Xm'][spotInd]
    y = spotData['Ym'][spotInd]
    frm = int(spotData['ome_idxs'][spotInd])
    eta, tth = sf.xyToEtaTthRecenter(x,y,params)
    
    # Track
    for t_ind, fnames in enumerate(dataFileSequence):
        if t_ind== 0:
            initSpot(spotInd,eta,tth,frm,fnames,trackFile,params)
        else:
            processSpot(spotInd,t_ind,fnames,trackFile,params)
            
def DoG(f,sigma,dsigma,gamma=2):
    g1 = gaussian_filter(f,sigma=sigma)
    g2 = gaussian_filter(f,sigma=sigma+dsigma)
    
    return (g2 - g1)/(sigma*dsigma) 
            
def detectBlobDoG(x):
    # 1. Compute normalized DoG apporximation of LoG
    dog_norm = DoG(x,sigma=2,dsigma=1.5)

    # 2. Pre-segmentation
    hess_mat = hessian_matrix(dog_norm)
    D1 = np.zeros(hess_mat[0].shape)
    D2 = np.zeros(hess_mat[0].shape)
    D3 = np.zeros(hess_mat[0].shape)
    for i1 in range(hess_mat[0].shape[0]):
        for i2 in range(hess_mat[0].shape[1]):
            for i3 in range(hess_mat[0].shape[2]):
                h_mat = np.array([[hess_mat[0][i1,i2,i3],hess_mat[1][i1,i2,i3],hess_mat[2][i1,i2,i3]],
                                  [hess_mat[1][i1,i2,i3],hess_mat[3][i1,i2,i3],hess_mat[4][i1,i2,i3]],
                                  [hess_mat[2][i1,i2,i3],hess_mat[4][i1,i2,i3],hess_mat[5][i1,i2,i3]]])
                D1[i1,i2,i3] = h_mat[0,0]
                D2[i1,i2,i3] = np.linalg.det(h_mat[:2,:2])
                D3[i1,i2,i3] = np.linalg.det(h_mat)

    negDefIndicator = (D1 > 0) & (D2 > 0) & (D3 > 0)
    blobs, num_blobs = label(negDefIndicator)
    return blobs, num_blobs, hess_mat

def blobFeaturesDoG(x,blobs,num_blobs,hess_mat):
    hess_T = []
    eigs_T = []
    RT_T = []
    ST_T = []
    AT_T = []
    for i in range(1,num_blobs+1):
        hess_i = []
        for component in hess_mat:
            hess_i.append(np.sum(component[blobs == i]))
        hess_T.append(hess_i)
        h_mat = np.array([[hess_i[0],hess_i[1],hess_i[2]],
                          [hess_i[1],hess_i[3],hess_i[4]],
                          [hess_i[2],hess_i[4],hess_i[5]]])
        eigs_i = np.linalg.eig(h_mat)
        eigs_T.append(eigs_i[0])
        
        # Blobness feature
        RT = 3*np.abs(eigs_i[0][0]*eigs_i[0][1]*eigs_i[0][2])**(2/3)/(np.abs(eigs_i[0][0]*eigs_i[0][1]) +
              np.abs(eigs_i[0][0]*eigs_i[0][2]) +
              np.abs(eigs_i[0][2]*eigs_i[0][1]))
        # Flatness feature
        ST = np.sqrt(eigs_i[0][0]**2 + eigs_i[0][1]**2 + eigs_i[0][2]**2)
        # Average intensity feature
        AT = np.mean(x[blobs == i])
        RT_T.append(RT)
        ST_T.append(ST)
        AT_T.append(AT)
    return RT_T,ST_T,AT_T
        
def gatherSimTruth(spotFiles,spotInds):
    spotDataList = []
    for i in range(len(spotFiles)):
        if os.path.exists(spotFiles[i]):
            spotDataList.append(np.load(spotFiles[i])) 
    numSpots = len(spotInds)
    spotIndArray = np.zeros((numSpots,5))
    spotIndArray[:,0] = np.arange(numSpots)
    
    for spotInd in range(numSpots):
        spot_id = sf.getSpotID(spotDataList[0],spotInd)
        for t in range(1,5):        
            m_ind = sf.matchSpotID(spotDataList[t],spot_id,spotInd)
            spotIndArray[spotInd,t] = m_ind
    return spotIndArray