#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spotfetch.visualization

Tools for visualizing X-ray diffraction data and results.

Created on: Fri Nov 15 23:00:28 2024
Author: Daniel Banco
"""
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import glob
import spotfetch as sf
import spotfetch.detectors as dt


def roiAdjacent(ind,tth,eta,frame,omFrms,timeFrms,params,dataDir):
    
    fig, axes = plt.subplots(len(omFrms), len(timeFrms), figsize=(len(omFrms), len(timeFrms)))
    
    # add enumerate HERR ND THEN coordinate subplots
    for i,frm in enumerate(omFrms):
        frm = sf.wrapFrame(frm)
        
        for j,t in enumerate(timeFrms):
            fnames = sf.timeToFile(t,dataDir)
            roi = dt.loadPolarROI(fnames,tth,eta,frm,params)
            img = axes[i,j].imshow(roi)
            fig.colorbar(img, ax=axes[i,j])
            axes[i,j].set_title(f'$\omega={frm}$, t={t}')

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
    ome_idxs = spotData['ome_idxs']-2
    # Get spot indices at frame
    if grains == []:
        spotInds = np.where(ome_idxs == frame)[0]
    else:
        spotInds = sf.findSpots(spotData,grains=grains,frm=frame)
    
    roiSize = params['roiSize']
    imSize = params['imSize']
    center = (imSize[0]//2,imSize[1]//2,)

    if os.path.isfile(fnames):
        fnames = [fnames]
    if detectFrame == []:
        detectFrame = frame
    if params['detector'] == 'dexela':
        b = dt.load_dex(fnames,params,detectFrame)
    elif params['detector'] == 'eiger':
        b = dt.load_eiger(fnames,params,detectFrame)
    elif params['detector'] == 'eiger_sim':
        b = dt.load_eiger_sim(fnames,params,detectFrame)
    
    
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
        detectDist, mmPerPixel, ff_trans = dt.loadYamlData(params,tth,eta)
        
        # 2. Construct rad, eta domain
        rad_dom, eta_dom = dt.polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
   
        rad = np.round(detectDist*np.tan(tth)/mmPerPixel)
        
        wedge = patches.Wedge([center[1],center[0]],rad+roiSize[0]/2,\
                180/np.pi*eta_dom[0],180/np.pi*eta_dom[-1],\
                linewidth=1,width=roiSize[1],fill=0,color='r')
            
        plt.gca().add_patch(wedge)
        # x = rad*np.cos(eta) + center[1]
        # y = rad*np.sin(eta) + center[0]
        # plt.plot(x,y,color='r',marker='x')
        
    return fig

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
                        omega = sf.frameToOmega(trackData[t][k][j]['frm'])
                        axes[k].scatter(t,omega,marker='s',color='b')
                        axes[k].set_title(f'spot {k}')
                        axes[k].set_xlim((0,20))

    return fig

def roiTrackVisual(spotInd,spotData,dome,num_cols,scanRange,dataPath,trackPath,truthPath,params):

    T = len(scanRange)
    N_ome = 2*dome+1
    
    num_figs = int(np.ceil(T/num_cols))
    
    fig_list = []
    axes_list = []
    
    # Create a figure and a set of subplots
    for i in range(num_figs):
        fig, axes = plt.subplots(N_ome, num_cols, figsize=(20, 15))
        
        # Remove x and y ticks for clarity
        for ax_row in axes:
            for ax in ax_row:
                ax.set_xticks([])
                ax.set_yticks([])
        
        fig_list.append(fig)
        axes_list.append(axes)
        
        i1 = 0 + i*num_cols
        i2 = num_cols-1 + i*num_cols
        j1 = scanRange[i1]
        if i2 >= T:
            j2 = scanRange[-1]
        else:
            j2 = scanRange[i2]
        # Add common labels
        fig_list[i].text(0.04, 0.5, r'$\omega$ frame', va='center', rotation='vertical', fontsize=24)
        fig_list[i].text(0.5, 0.04, 'Scan #', ha='center', fontsize=24)
        fig_list[i].text(0.5, 0.95, f'Spot {spotInd}, Scans {j1}-{j2}', ha='center', fontsize=32)      
        fig_list[i].subplots_adjust(wspace=0.05, hspace=0.01)
    
    # Load track data
    track_file = os.path.join(trackPath,f'trackData_{spotInd}.pkl')
    if os.path.exists(track_file):
        with open(track_file, 'rb') as f:
            track_data = pickle.load(f)
            print('Track loaded')
    else:
        track_data = []
        trackFound = False
        
    # Load truth data 
    truth_file = os.path.join(trackPath,f'truthData_{spotInd}.pkl')
    if os.path.exists(truth_file):
        with open(truth_file, 'rb') as f:
            truth_data = pickle.load(f)
            print('Truth loaded')
    else:   
        T = len(scanRange)
        truth_data = []
        truthFound = False
    
    eta0 = spotData['etas'][spotInd]
    tth0 = spotData['tths'][spotInd]
    frm0 = spotData['ome_idxs'][spotInd]
    etaRoi = eta0
    tthRoi = tth0
    frmRange = np.arange(frm0-dome,frm0+dome+1)
    
    for scan_ind in range(T):
        for om_ind in range(N_ome):
            i = int(np.floor(scan_ind/num_cols))
            j = np.mod(scan_ind,num_cols)
            ax = axes_list[i][om_ind,j]
            scan = scanRange[scan_ind]
            frm = sf.wrapFrame(frmRange[om_ind])
            print(f'Scan {scan}, Frame {frm}')
            # 1. Show ROI
            sf.showROI(ax,dataPath, scan, frm,\
                    tthRoi, etaRoi, params)
            sf.showInitial(ax,etaRoi,tthRoi,eta0,tth0,params)
            # 2. Show the track and truth if exist
            if len(track_data) > 0:
                [trackFound, om_ind1] = sf.checkTrack(track_data,scan_ind,scan,frm)
            if len(truth_data) > 0:
                [truthFound, om_ind2] = sf.checkTruth(truth_data,scan_ind,scan,frm)      
            
            # Add tracks and truth if available
            if trackFound:
                track = track_data[scan_ind][om_ind1].copy()
                sf.showTrack(ax,track,etaRoi,tthRoi,\
                          eta0,tth0,params)
            if truthFound:
                truth = truth_data[scan_ind][om_ind2].copy()
                sf.showTruth(ax,truth,etaRoi,tthRoi,params)
            
            # Label plots
            if om_ind == N_ome-1:
                ax.set_xlabel(f'{scan}')
            if scan_ind == 0:
                ax.set_ylabel(f'{frm}')


def roiTrackVisual2(spotInds,spotData,dome,scanRange,trackPath,dataPath,params):
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
            omRange[i] = sf.wrapFrame(om)
        
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
        if len(track_data[0]) == 0:
            continue
        track = track_data[0]
        eta0 = track[0]['etaRoi']
        tth0 = track[0]['tthRoi']
        
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
                fnames = sf.timeToFile(scanRange[j],dataPath)
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
                    template = dataPath.format(num2=scan)
                    fnames = glob.glob(template)
                    roi = dt.loadPolarROI(fnames,tthRoi,etaRoi,frm,params)
                else:
                    roi = dt.loadPolarROI(fnames,tthRoi,etaRoi,frm,params)
                ax.imshow(roi)
                ax.text(1, 4, f'{roi.max():.2f}', color='white', fontsize=12, weight='bold')
                
                if (p1 > 0) & (p2 > 0):
                    # Plot rect
                    eta = etaTrack[i,j]
                    tth = tthTrack[i,j]
                    detectDist, mmPerPixel, ff_trans = dt.loadYamlData(params,tth,eta)
                    boxSize = 10
                    start_col = round(p1 - boxSize//2)
                    start_row = round(p2 - boxSize//2)
                    rect = plt.Rectangle((start_col, start_row), boxSize, boxSize,
                                          linewidth=1, edgecolor='r', facecolor='none')
                    ax.add_patch(rect)
                    
                    # Show Previous Track
                    rad_dom, eta_dom = dt.polarDomain(detectDist, mmPerPixel,\
                                          tthRoi, etaRoi, params['roiSize'])
                    radRoi = detectDist*np.tan(tthRoi)/mmPerPixel
                    
                    y_pos = (radRoi-rad_dom[0])/(rad_dom[-1]-rad_dom[0])*39
                    x_pos = (etaRoi-eta_dom[0])/(eta_dom[-1]-eta_dom[0])*39
                    ax.plot(x_pos,y_pos,marker='o',markersize=16,\
                            fillstyle='none',color='red')
                    
                    rad0 = detectDist*np.tan(tth0)/mmPerPixel
                    y_pos = (rad0-rad_dom[0])/(rad_dom[-1]-rad_dom[0])*39
                    x_pos = (eta0-eta_dom[0])/(eta_dom[-1]-eta_dom[0])*39
                    ax.plot(x_pos,y_pos,marker='x',markersize=10,color='yellow')
                    
                    
                if i == rows-1:
                    ax.set_xlabel(f'{scan}')
                if j == 0:
                    ax.set_ylabel(f'{frm}')
                    
                ax.set_xticks([])
                ax.set_yticks([])
        plt.draw()
        # plt.savefig(f'fig_{k}.png')
        # plt.close(fig)