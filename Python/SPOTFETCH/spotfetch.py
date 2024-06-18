# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 13:22:47 2024

@author: dpqb1
"""
import numpy as np
import scipy as sp
import h5py
# import hdf5plugin
import pickle
import pandas as pd
import os
import subprocess
import time
import yaml
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from hexrd.fitting import fitpeak
# from hexrd.fitting import peakfunctions as pkfuncs

def omegaToFrame(omega,startFrame=4,endFrame=1441,omegaRange=360,startOmeg = 0):
    step = omegaRange/endFrame
    frame = (np.floor((omega-startOmeg)/step) + startFrame).astype(int)
    return frame

def frameToOmega(frame,startFrame=4,endFrame=1441,omegaRange=360,startOmeg = 0):
    step = omegaRange/endFrame
    omega = (frame - startFrame)*step + startOmeg
    return omega

def read_yaml(file_path):
    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)
    return data

def loadYamlData(params,tth,eta):
    yamlFile = params['yamlFile']
    if params['detector'] == 'dexela':
        return loadYamlDataDexela(yamlFile,tth,eta)
    elif params['detector'] == 'eiger':
        return loadYamlDataEiger(yamlFile,tth,eta)

def loadYamlDataDexela(yamlFile,tth,eta):
    yamlData = read_yaml(yamlFile)
    ff1_trans = yamlData['detectors']['ff1']['transform']['translation']
    ff2_trans = yamlData['detectors']['ff2']['transform']['translation']
    ff_trans = [ff1_trans[0],ff1_trans[1],ff2_trans[0],ff2_trans[1]]
    if abs(eta)>(np.pi/2):
        detectDist = -ff2_trans[2]
        mmPerPixel = yamlData['detectors']['ff2']['pixels']['size'][0]
    else:
        detectDist = -ff1_trans[2] 
        mmPerPixel = yamlData['detectors']['ff1']['pixels']['size'][0]
    
    return detectDist, mmPerPixel,ff_trans

def loadYamlDataEiger(yamlFile,tth,eta):
    yamlData = read_yaml(yamlFile)
    ff_trans = yamlData['detectors']['eiger']['transform']['translation']
    detectDist = ff_trans[2]
    mmPerPixel = yamlData['detectors']['eiger']['pixels']['size'][0]

    return detectDist, mmPerPixel, ff_trans

def load_eiger(fname,yamlFile,t,t2=None,imSize=(4362,4148)):

    if t2 is None:
        with h5py.File(fname, 'r') as file1:
            group = file1['entry']
            data1 = group['data']
            data2 = data1['data']
            img = data2[t,:,:]
    else:
        with h5py.File(fname, 'r') as file1:
            group = file1['entry']
            data1 = group['data']
            data2 = data1['data']
            img = data2[t:t2,:,:]

        # Take the maximum of multiple frames in time
        img = np.max(img, axis=0)


    # Maybe just need to adjust so center lies at center of image
    
    return img    
    
def load_dex(fname1,fname2,params,t,t2=None,dexSize=(3888,3072)):
    if t2 is None:
        with h5py.File(fname1, 'r') as file1, h5py.File(fname2, 'r') as file2:
            img1 = file1['/imageseries/images'][t, :, :]
            img2 = file2['/imageseries/images'][t, :, :]
    else:
        with h5py.File(fname1, 'r') as file1, h5py.File(fname2, 'r') as file2:
            img1 = file1['/imageseries/images'][t:t2, :, :]
            img2 = file2['/imageseries/images'][t:t2, :, :]

        # Take the maximum of multiple frames in time
        img1 = np.max(img1, axis=0)
        img2 = np.max(img2, axis=0)

    imSize = params['imSize']
    yamlFile = params['yamlFile']
    # Pad image
    bpad = np.zeros(imSize)
    center = (imSize[0]/2, imSize[1]/2)
    
    # Shift each panel
    yamlData = read_yaml(yamlFile)
    mmPerPixel  = yamlData['detectors']['ff2']['pixels']['size'][0]
    ff1_trans = yamlData['detectors']['ff1']['transform']['translation']
    ff2_trans = yamlData['detectors']['ff2']['transform']['translation']
    ff1_tx = ff1_trans[0]
    ff1_ty = ff1_trans[1]
    ff2_tx = ff2_trans[0]
    ff2_ty = ff2_trans[1]
    ff2_xshift = int(round(ff2_tx/mmPerPixel))
    ff1_xshift = int(round(ff1_tx/mmPerPixel))
    ff2_yshift = int(round(ff2_ty/mmPerPixel))
    ff1_yshift = int(round(ff1_ty/mmPerPixel))
    
    # negative sign on y shift because rows increase downwards
    ff1r1 = int(center[0]-dexSize[0]/2-ff1_yshift)
    ff1r2 = int(center[0]+dexSize[0]/2-ff1_yshift)
    ff2r1 = int(center[0]-dexSize[0]/2-ff2_yshift)
    ff2r2 = int(center[0]+dexSize[0]/2-ff2_yshift)
    
    ff1c1 = int(center[1]-dexSize[1]/2+ff1_xshift)
    ff1c2 = int(center[1]+dexSize[1]/2+ff1_xshift)
    ff2c1 = int(center[1]-dexSize[1]/2+ff2_xshift)
    ff2c2 = int(center[1]+dexSize[1]/2+ff2_xshift)
        
    bpad[ff2r1:ff2r2,ff2c1:ff2c2] = np.flipud(img2)
    bpad[ff1r1:ff1r2,ff1c1:ff1c2] = np.fliplr(img1)
    
    return bpad
def getInterpParams(tth,eta,params):
    if params['detector'] == 'dexela':
        return getInterpParamsDexela(tth,eta,params)
    elif params['detector'] == 'eiger':
        return getInterpParamsEiger(tth,eta,params)

def getInterpParamsDexela(tth,eta,params):
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    imSize = params['imSize']
    
    center = (imSize[0]/2,imSize[1]/2)
    detectDist, mmPerPixel, ff_trans = loadYamlDataDexela(yamlFile,tth,eta)
    
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    x_cart, y_cart = fetchCartesian(rad_dom,eta_dom,center)
    ff1_pix, ff2_pix = panelPixels(ff_trans,mmPerPixel,imSize)
    
    new_center = np.array([center[0] - y_cart[0], center[1] - x_cart[0]])
    roiShape = getROIshape(x_cart, y_cart, ff1_pix, ff2_pix, center)
    
    Ainterp = bilinearInterpMatrix(roiShape,rad_dom,eta_dom,new_center)
    
    return Ainterp, new_center, x_cart, y_cart

def getInterpParamsEiger(tth,eta,params):
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    imSize = params['imSize']
    
    center = (imSize[0]/2,imSize[1]/2)
    detectDist, mmPerPixel, ff_trans = loadYamlDataEiger(yamlFile,tth,eta)
    
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    x_cart, y_cart = fetchCartesian(rad_dom,eta_dom,center)
    ff1_pix, ff2_pix = panelPixels(ff_trans,mmPerPixel,imSize)
    
    new_center = np.array([center[0] - y_cart[0], center[1] - x_cart[0]])
    roiShape = getROIshape(x_cart, y_cart, ff1_pix, ff2_pix, center)
    
    Ainterp = bilinearInterpMatrix(roiShape,rad_dom,eta_dom,new_center)
    
    return Ainterp, new_center, x_cart, y_cart

def loadPolarROI(fname1,fname2,tth,eta,frame,params):
    if params['detector'] == 'dexela':
        roi = loadDexPolarRoi(fname1,fname2,tth,eta,frame,params)
    elif params['detector'] == 'eiger':
        roi = loadEigerPolarRoi(fname1,fname2,tth,eta,frame,params)
    return roi
    
def loadDexPolarRoi(fname1,fname2,tth,eta,frame,params):
    # 0. Load params, YAML data
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans = loadYamlDataDexela(yamlFile,tth,eta)
    
    # 1. Construct rad, eta domain
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    
    # 2. Construct interpolation matrix
    Ainterp,new_center,x_cart,y_cart = getInterpParamsDexela(tth,eta,params)
    
    # 3. Load needed Cartesian ROI pixels
    ff1_pix, ff2_pix = panelPixels(ff_trans,mmPerPixel)
    roi = loadDexPanelROI(x_cart,y_cart,ff1_pix,ff2_pix,fname1,fname2,frame,params)
    
    # 4. Apply interpolation matrix to Cartesian pixels get Polar values
    roi_polar_vec = Ainterp.dot(roi.flatten())
    
    # 5. Reshape and output roi
    roi_polar = np.reshape(roi_polar_vec,(len(rad_dom),len(eta_dom)))
    
    return roi_polar

def loadEigerPolarRoi(fname1,fname2,tth,eta,frame,params):
    # 0. Load params, YAML data
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    imSize = params['imSize']
    detectDist, mmPerPixel, ff_trans = loadYamlData(yamlFile,tth,eta)
    
    # 1. Construct rad, eta domain
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    
    # 2. Construct interpolation matrix
    Ainterp,new_center,x_cart,y_cart = getInterpParamsDexela(imSize,tth,eta,roiSize,yamlFile)
    
    # 3. Load needed Cartesian ROI pixels
    ff1_pix, ff2_pix = panelPixels(ff_trans,mmPerPixel)
    roi = loadDexPanelROI(x_cart,y_cart,ff1_pix,ff2_pix,fname1,fname2,frame,imSize)
    
    # 4. Apply interpolation matrix to Cartesian pixels get Polar values
    roi_polar_vec = Ainterp.dot(roi.flatten())
    
    # 5. Reshape and output roi
    roi_polar = np.reshape(roi_polar_vec,(len(rad_dom),len(eta_dom)))
    
    return roi_polar

def loadDexPolarRoiPrecomp(fname1,fname2,yamlFile,tth,eta,frame,\
                    imSize,roiSize,interp_params):
    # 0. Load YAML data
    detectDist, mmPerPixel, ff_trans = loadYamlData(yamlFile,tth,eta)
    
    # 1. Construct rad, eta domain
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    
    # 2. Load precomputed interpolation data
    Ainterp = interp_params['A']
    x_cart = interp_params['x_cart']
    y_cart = interp_params['y_cart']
    
    # 3. Load needed Cartesian ROI pixels
    ff1_pix, ff2_pix = panelPixels(ff_trans,mmPerPixel)
    roi = loadDexPanelROI(x_cart,y_cart,ff1_pix,ff2_pix,fname1,fname2,frame,imSize)
    
    # 4. Apply interpolation matrix to Cartesian pixels get Polar values
    roi_polar_vec = Ainterp.dot(roi.flatten())
    
    # 5. Reshape and output roi
    roi_polar = np.reshape(roi_polar_vec,(len(rad_dom),len(eta_dom)))
    
    return roi_polar
          
def loadSpotsAtFrame(spot_data,fname1,fname2,frame,params):
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
        roi_polar = loadDexPolarRoi(fname1,fname2,tth,eta,frame,params)
        # roi_polar_p1 = loadDexPolarRoi(fname1,fname2,yamlFile,\
        #             tth,eta,frame+1,imSize,roi_size)
        # roi_polar_m1 = loadDexPolarRoi(fname1,fname2,yamlFile,\
        #             tth,eta,frame-1,imSize,roi_size)
            
        roi_list.append(roi_polar)
        
    return roi_list

def loadDexPanelROI(x_cart,y_cart,ff1_pix,ff2_pix,fname1,fname2,frame,params,\
                    dexShape=(3888,3072)):
    imSize = params['imSize']
    center = (imSize[0]/2,imSize[1]/2)
    
    if x_cart[0] < ff2_pix[0]: x_cart[0] = ff2_pix[0]
    if x_cart[1] > ff1_pix[1]: x_cart[1] = ff1_pix[1]
    # Determine which panel we are on
    if x_cart[0] < center[1]:
        if y_cart[0] < ff2_pix[2]: y_cart[0] = ff2_pix[2]
        if y_cart[1] > ff2_pix[3]: y_cart[1] = ff2_pix[3]
        x_pan = x_cart - ff2_pix[0]
        y_pan = y_cart - ff2_pix[2]
        # Need to account for flip up down
        midLine = (dexShape[0]-1)/2
        flip0 = y_pan[0] + 2*(midLine-y_pan[0])
        flip1 = y_pan[1] + 2*(midLine-y_pan[1])
        y_pan[0] = min(flip0,flip1)
        y_pan[1] = max(flip0,flip1)
        with h5py.File(fname2, 'r') as file:
            img = file['/imageseries/images'][frame, y_pan[0]:y_pan[1],x_pan[0]:x_pan[1]]
            img = np.flipud(img)
    elif x_cart[0] > center[1]:
        if y_cart[0] < ff1_pix[2]: y_cart[0] = ff1_pix[2]      
        if y_cart[1] > ff1_pix[3]: y_cart[1] = ff1_pix[3]
        x_pan = x_cart - ff1_pix[0]
        y_pan = y_cart - ff1_pix[2]
        # Need to account for flip left right
        midLine = (dexShape[1]-1)/2
        flip0 = x_pan[0] + 2*(midLine-x_pan[0])
        flip1 = x_pan[1] + 2*(midLine-x_pan[1])
        x_pan[0] = min(flip0,flip1)
        x_pan[1] = max(flip0,flip1)
        with h5py.File(fname1, 'r') as file:
            img = file['/imageseries/images'][frame, y_pan[0]:y_pan[1],x_pan[0]:x_pan[1]]
            img = np.fliplr(img)
      
    return img

def getROIshape(x_cart,y_cart,ff1_pix,ff2_pix,center,dexShape=(3888,3072)):
    
    if x_cart[0] < ff2_pix[0]: x_cart[0] = ff2_pix[0]
    if x_cart[1] > ff1_pix[1]: x_cart[1] = ff1_pix[1]
    # Determine which panel we are on
    if x_cart[0] < center[1]:
        if y_cart[0] < ff2_pix[2]: y_cart[0] = ff2_pix[2]
        if y_cart[1] > ff2_pix[3]: y_cart[1] = ff2_pix[3]
        x_pan = x_cart - ff2_pix[0]
        y_pan = y_cart - ff2_pix[2]
        # Need to account for flip up down
        midLine = (dexShape[0]-1)/2
        flip0 = y_pan[0] + 2*(midLine-y_pan[0])
        flip1 = y_pan[1] + 2*(midLine-y_pan[1])
        y_pan[0] = min(flip0,flip1)
        y_pan[1] = max(flip0,flip1)
    elif x_cart[0] > center[1]:
        if y_cart[0] < ff1_pix[2]: y_cart[0] = ff1_pix[2]      
        if y_cart[1] > ff1_pix[3]: y_cart[1] = ff1_pix[3]
        x_pan = x_cart - ff1_pix[0]
        y_pan = y_cart - ff1_pix[2]
        # Need to account for flip left right
        midLine = (dexShape[1]-1)/2
        flip0 = x_pan[0] + 2*(midLine-x_pan[0])
        flip1 = x_pan[1] + 2*(midLine-x_pan[1])
        x_pan[0] = min(flip0,flip1)
        x_pan[1] = max(flip0,flip1)
    
    return (y_pan[1]-y_pan[0],x_pan[1]-x_pan[0])

def bilinearInterpMatrix(roiShape,rad_dom,eta_dom,center):
    Ainterp = np.zeros((len(rad_dom)*len(eta_dom),roiShape[0]*roiShape[1]))
    k = 0
    for r in rad_dom:
        for eta in eta_dom:
            x = r*np.cos(eta) + center[1]
            y = -r*np.sin(eta) + center[0]
            x1 = np.floor( x );
            x2 = np.ceil(  x );
            y1 = np.floor( y );
            y2 = np.ceil(  y );
            
            Ainterp[k,int(x2 + y2*roiShape[1])] = (x-x1)*(y-y1)
            Ainterp[k,int(x2 + y1*roiShape[1])] = (x2-x)*(y-y1)
            Ainterp[k,int(x1 + y2*roiShape[1])] = (x-x1)*(y2-y)
            Ainterp[k,int(x1 + y1*roiShape[1])] = (x2-x)*(y2-y)

            k += 1
    Ainterp = sp.sparse.coo_array(Ainterp)
    return Ainterp

def panelPixels(ff_trans,mmPerPixel,imSize=(4888,7300),dexShape=(3888,3072)):
    center = (imSize[0]/2,imSize[1]/2)
    ff1_tx = ff_trans[0]
    ff1_ty = ff_trans[1]
    ff2_tx = ff_trans[2]
    ff2_ty = ff_trans[3]
    ff2_xshift = int(round(ff2_tx/mmPerPixel))
    ff1_xshift = int(round(ff1_tx/mmPerPixel))
    ff2_yshift = int(round(ff2_ty/mmPerPixel))
    ff1_yshift = int(round(ff1_ty/mmPerPixel))
    
    # negative sign on y shift because rows increase downwards
    ff1r1 = int(center[0]-dexShape[0]/2-ff1_yshift)
    ff1r2 = int(center[0]+dexShape[0]/2-ff1_yshift)
    ff2r1 = int(center[0]-dexShape[0]/2-ff2_yshift)
    ff2r2 = int(center[0]+dexShape[0]/2-ff2_yshift)
    
    ff1c1 = int(center[1]-dexShape[1]/2+ff1_xshift)
    ff1c2 = int(center[1]+dexShape[1]/2+ff1_xshift)
    ff2c1 = int(center[1]-dexShape[1]/2+ff2_xshift)
    ff2c2 = int(center[1]+dexShape[1]/2+ff2_xshift)
    
    ff1_pixels = [ff1c1,ff1c2,ff1r1,ff1r2]
    ff2_pixels = [ff2c1,ff2c2,ff2r1,ff2r2]
    return ff1_pixels, ff2_pixels

def polarDomain(detectDist,mmPerPixel,tth,eta,roi_size):
    r = np.round(detectDist*np.tan(tth)/mmPerPixel)
    r1 = r - roi_size//2
    r2 = r + roi_size//2

    rad_domain = np.arange(r1,r2,1)
    
    deta = 1/r2
    eta1 = eta - roi_size/2*deta
    eta2 = eta + roi_size/2*deta
    eta_domain = np.linspace(eta1,eta2,roi_size)
    return rad_domain, eta_domain

def fetchCartesian(rad_dom,eta_dom,center):
    rad1 = rad_dom[0]*np.cos(eta_dom)
    rad2 = rad_dom[-1]*np.cos(eta_dom)
    x_min1 = (np.min(rad1 + center[1]))
    x_max1 = (np.max(rad2 + center[1]))
    x_max2 = (np.max(rad1 + center[1]))    
    x_min2 = (np.min(rad2 + center[1]))
    x_min = int(np.floor(min(x_min1,x_min2)))
    x_max = int(np.ceil(max(x_max1,x_max2)))
    
    rad1 = rad_dom[0]*np.sin(eta_dom)
    rad2 = rad_dom[-1]*np.sin(eta_dom)
    y_min1 = (np.min(-rad1 + center[0]))
    y_max1 = (np.max(-rad2 + center[0]))
    y_max2 = (np.max(-rad1 + center[0]))
    y_min2 = (np.min(-rad2 + center[0]))
    y_min = int(np.floor(min(y_min1,y_min2)))
    y_max = int(np.ceil( max(y_max1,y_max2)))
    
    return np.array([x_min,x_max+1]),np.array([y_min,y_max+1])


def plot_ring(radius,center):
    eta = np.linspace(0,2*np.pi,600)
    x = radius*np.sin(eta) + center[1]
    y = radius*np.cos(eta) + center[0]
    plt.plot(x, y, color='red', linewidth=1)

def add_spot_rect(params,tth,eta,roi_size=40):
    detectDist,mmPerPixel,ff_trans = loadYamlData(params, tth, eta) 
    imSize = params['imSize']
    center = (imSize[0]//2,imSize[1]//2)
    radius = detectDist*np.tan(tth)/mmPerPixel
    x = radius*np.cos(eta) + center[1]
    y = -radius*np.sin(eta) + center[0] #Negative sign because rows increase downwards?
    
    # Calculate the starting and ending indices for rows and columns
    start_col = round(x - roi_size//2)
    start_row = round(y - roi_size//2)
    
    rect = plt.Rectangle((start_col, start_row), roi_size, roi_size,
                              linewidth=1, edgecolor='r', facecolor='none')
    plt.gca().add_patch(rect)
    
def get_roi(b,detectDist,mmPerPixel,center,tth,eta,roi_size=40):
    radius = detectDist*np.tan(tth)/mmPerPixel
    x = radius*np.cos(eta) + center[1]
    y = -radius*np.sin(eta) + center[0] #Negative sign because rows increase downwards?
    
    # Calculate the starting and ending indices for rows and columns
    start_col = round(x - roi_size//2)
    end_col = round(x + roi_size//2)
    start_row = round(y - roi_size//2)
    end_row = round(y + roi_size//2)
    
    roi = b[start_row:end_row,start_col:end_col]
    return roi

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
    
def plotSpotWedges(spotData,fname1,fname2,frame,params):
    tths = spotData['tths']
    etas = spotData['etas']     
    ome_idxs = spotData['ome_idxs'] 
    # Get spot indices at frame
    spotInds = np.where(ome_idxs == frame)[0]
    
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    imSize = params['imSize']
    center = (imSize[0]//2,imSize[1]//2,)

    b = load_dex(fname1,fname2,params,frame)
    
    fig = plt.figure()
    plt.imshow(b)
    plt.clim(290,550)
    
    for ind in spotInds:
        # 0. Load spot information
        tth = tths[ind]
        eta = -etas[ind]

        # 1. Load YAML data
        detectDist, mmPerPixel, ff_trans = loadYamlData(params,tth,eta)
        
        # 2. Construct rad, eta domain
        rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
   
        rad = np.round(detectDist*np.tan(tth)/mmPerPixel)
        
        wedge = patches.Wedge([center[1],center[0]],rad+roiSize/2,\
                180/np.pi*eta_dom[0],180/np.pi*eta_dom[-1],\
                linewidth=1,width=roiSize,fill=0,color='r')
        plt.gca().add_patch(wedge)
        x = rad*np.cos(eta) + center[1]
        y = rad*np.sin(eta) + center[0]
        plt.plot(x,y,color='r',marker='x')
        
    return fig

def plotSpotRectangles(spot_data,fname1,fname2,\
                   yamlFile,frame,center,roi_size):
    tths = spot_data['tths']
    etas = spot_data['etas']     
    ome_idxs = spot_data['ome_idxs'] 
    spotInds = np.where(ome_idxs == frame)[0]
    
    b = load_dex(fname1,fname2,yamlFile,frame)
    
    fig = plt.figure()
    plt.imshow(b)
    plt.clim(290,550)
        
    for ind in spotInds:
        # 0. Load spot information
        tth = tths[ind]
        eta = etas[ind]
    
        # 1. Load YAML data
        detectDist, mmPerPixel, ff_trans = loadYamlData(yamlFile,tth,eta)
        
        # 2. Add rectangle
        add_spot_rect(detectDist,mmPerPixel,center,tth,eta)       
        
        rad = np.round(detectDist*np.tan(tth)/mmPerPixel)
        x = rad*np.cos(eta) + center[1]
        y = -rad*np.sin(eta) + center[0]
        
        plt.plot(x,y,color='r',marker='x')
        
    return fig

def wrapFrame(frm,frm0=4,frmEnd=1441):
    return np.mod(frm-frm0,frmEnd)+frm0

def timeToFile(t,fDir):
    dNum = str(t)
    topDir = os.path.join(fDir + dNum, 'ff')
    fnames = os.listdir(topDir)
    fnames.sort()
    fname1 = os.path.join(topDir,fnames[0])
    fname2 = os.path.join(topDir,fnames[1])
    
    return fname1, fname2

def pathToFile(path):
    fnames = os.listdir(path)
    fnames.sort()
    return fnames

def timeToFileB(t):
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
        
    return fname1,fname2

def roiAdjacent(yamlFile,ind,tth,eta,imSize,roi_size,frame,omFrms,timeFrms):
    
    fig, axes = plt.subplots(len(omFrms), len(timeFrms), figsize=(len(omFrms), len(timeFrms)))
    
    # add enumerate HERR ND THEN coordinate subplots
    for i,frm in enumerate(omFrms):
        frm = wrapFrame(frm)
        
        for j,t in enumerate(timeFrms):
            fname1,fname2 = timeToFileB(t)
            roi = loadDexPolarRoi(fname1,fname2,yamlFile,tth,eta,frm,imSize,roi_size)
            img = axes[i,j].imshow(roi)
            fig.colorbar(img, ax=axes[i,j])
            axes[i,j].set_title(f'$\omega={frm}$, t={t}')
            
    
def evaluateROI(fname1,fname2,prevTracks,tth,eta,frm,scan,params):
    # 0. Parameters
    roiSize = params['roiSize']
    
    # 1. Load ROI
    roi = loadPolarROI(fname1,fname2,tth,eta,frm,params)
    tth_vals, eta_vals = np.indices(roi.shape)
    
    # 2. Estimate peak parameters (use from previous timestep)
    p0 = fitpeak.estimate_pk_parms_2d(eta_vals,tth_vals,roi,"gaussian")
    
    # 3. Fit peak
    p = fitpeak.fit_pk_parms_2d(p0,eta_vals,tth_vals,roi,"gaussian")
    
    if (p[1] > roiSize-0.5) | (p[2] > roiSize-0.5) | (p[1] < -0.5) | (p[2] < -0.5):
        peakFound = False
        # print('Mean exceeds ROI')
        return 0, peakFound
    
    residual = fitpeak.fit_pk_obj_2d(p,eta_vals,tth_vals,roi,"gaussian")
    rel_error = np.linalg.norm(residual)/np.linalg.norm(roi.flatten())
    
    detectDist, mmPerPixel, ff_trans = loadYamlData(params,tth,eta)
    rad_dom, eta_dom = polarDomain(detectDist,mmPerPixel,tth,eta,roiSize)  
    
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
    detectDist, mmPerPixel, ff_trans = loadYamlData(params,tth,eta)
    rad_dom, eta_dom = polarDomain(detectDist,mmPerPixel,tth,eta,roiSize)  
    deta = eta_dom[1]-eta_dom[0]
    dtth = np.arctan(mmPerPixel/detectDist)
    
    if (p[1] > roiSize) | (p[2] > roiSize):
        peakFound = False
        # print('Mean exceeds ROI')
        return peakFound
    
    if (p[3]==0): p[3] += 0.001
    if (p[4]==0): p[4] += 0.001
    
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

def collectSpotsData(topPath,spotsDir):
    # Get a list of file names in the directory and sort them
    spotsPath = topPath + spotsDir
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
                
    invalid_nums = id_nums == -999        
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
    
    np.savez(topPath + spotsDir + '.npz',Xs=Xs,Ys=Ys,id_nums=id_nums,\
    tths=tths,etas=etas,omes=omes,ome_idxs=ome_idxs,grain_nums=grain_nums)

def spotTracker(dataPath,outPath,exsituPath,spotData,spotInds,params,scan1):
    # Initialize 
    initData = {}
    initData['tths'] = spotData['tths'][spotInds]
    initData['etas'] = spotData['etas'][spotInds]
    initData['frms'] = spotData['ome_idxs'][spotInds]
    
    # Initialize tracks
    fnames = pathToFile(exsituPath)
    fname1 = exsituPath+fnames[0]
    fname2 = exsituPath+fnames[1]
    i = 0
    t = 0
     # Scan index
    print('')
    print(f'Scan {t}, Spot:', end=" ")
    for k,s in enumerate(spotInds):
        print(f'{k}', end=" ")
        etaRoi = initData['etas'][s]
        tthRoi = initData['tths'][s]
        frm = initData['frms'][s]
        initSpot(k,s,t,etaRoi,tthRoi,frm,params,outPath,fname1,fname2)
    
    i += 1
    t = scan1
    while True:
        # Try reading in file for new scan
        dataDir = os.path.join(dataPath,f'{t}','ff')
        try:
            fnames = pathToFile(dataDir)
        except:
            print('No new data')
            time.sleep(1)
            continue
        # if read is successful...
        if len(fnames) > 0:
            fname1 = os.path.join(dataDir,fnames[0])
            fname2 = os.path.join(dataDir,fnames[1])

            print('')
            print(f'Scan {t}, Spot:', end=" ")
            for k,s in enumerate(spotInds): 
                print(f'{k}', end=" ")
                processSpot(k,s,t,params,outPath,fname1,fname2)
            
            i += 1
            t += 1
        
def initSpot(k,s,t,etaRoi,tthRoi,frm,params,outPath,fname1,fname2):
    prevTracks = []
    newTrack, peakFound = evaluateROI(fname1,fname2,prevTracks,\
                        tthRoi,etaRoi,int(frm),t,params)
    trackData = []
    trackData.append([newTrack])
    outFile = os.path.join(outPath,'outputs',f'trackData_{k}.pkl')
    with open(outFile, 'wb') as f:
        pickle.dump(trackData, f)
    
def processSpot(k,s,t,params,outPath,fname1,fname2):
    outFile = os.path.join(outPath,'outputs',f'trackData_{k}.pkl')
    # Load in track so far
    with open(outFile, 'rb') as f:
        trackData = pickle.load(f)
    trackData.append([])
    T = len(trackData)

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

        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fname1,fname2,prevTracks,\
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
    if len(trackData[T-1]) == 0:
        frm1 = prevTracks[0]['frm']
        frm2 = prevTracks[-1]['frm']
        expandRange = list(range(frm1-3,frm1)) + list(range(frm2+1,frm2+4))
        for frm in expandRange:
            frm = int(wrapFrame(frm))
            newTrack, peakFound = evaluateROI(fname1,fname2,prevTracks,\
                                tth,eta,frm,t,params)
            if peakFound: 
                # print(f'Peak found at frame {frm}')
                trackData[T-1].append(newTrack)
                compareTrack = trackData[T-1]
                frm1 = trackData[T-1][0]['frm']
                frm2 = trackData[T-1][-1]['frm']
                break
    
    # Incremental Search if we have a peak found
    # Search down
    if len(trackData[T-1]) > 0: peakFound = True
    while peakFound:
        frm1 = frm1 - 1
        frm = int(wrapFrame(frm1))
        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fname1,fname2,compareTrack,\
                                          tth,eta,frm,t,params)
        # Add to list if peakFound
        if peakFound: 
            # print(f'Found more at {frm1}')
            trackData[T-1].insert(0,newTrack)

    # Search up
    if len(trackData[T-1]) > 0: peakFound = True
    while peakFound:
        frm2 = frm2 + 1
        frm = int(wrapFrame(frm2))
        # Load ROI and fit peak
        newTrack, peakFound = evaluateROI(fname1,fname2,compareTrack,\
                                          tth,eta,frm,t,params)
        # Add to list if peakFound
        if peakFound: 
            # print(f'Found more at {frm2}')
            trackData[T-1].append(newTrack)

    with open(outFile, 'wb') as f:
        pickle.dump(trackData, f)
        
def processSpotJob(inputFile):

    with open(inputFile, 'rb') as f:
        k,s,t,params,outPath,fname1,fname2 = pickle.load(f)
        
    processSpot(k,s,t,params,outPath,fname1,fname2)
    
def collect_job_ids():
    """
    Collects job IDs from the output of the qstat command.

    Returns:
        list: A list of job IDs currently in the queue.
    """
    try:
        # Run the qstat command
        result = subprocess.run(['qstat'], capture_output=True, text=True, check=True)
        
        # Split the output into lines
        lines = result.stdout.strip().split('\n')
        
        # Extract job IDs from the lines (assuming the job ID is in the first column)
        job_ids = []
        for line in lines[2:]:  # Skip the header lines
            parts = line.split()
            if len(parts) > 0:
                job_id = parts[0]
                job_ids.append(job_id)
        
        return job_ids
    
    except subprocess.CalledProcessError as e:
        print(f"Error running qstat: {e}")
        return []
    
def check_job_status(job_id):
    result = subprocess.run(['qstat'], capture_output=True, text=True)
    return job_id in result.stdout

def wait_for_jobs(job_ids):
    while True:
        all_completed = True
        for job_id in job_ids:
            if check_job_status(job_id):
                all_completed = False
                break
        if all_completed:
            break
        time.sleep(5)  # Wait for 1 minute before checking again
    
        
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
    fname1 = os.path.join(exsituPath,fnames[0])
    fname2 = os.path.join(exsituPath,fnames[1])
    i = 0
    t = scan1  # Scan index
    print('')
    print(f'Scan {t}, Spot:', end=" ")
    for k, s in enumerate(spotInds):
        print(f'{k}', end=" ")
        etaRoi = initData['etas'][s]
        tthRoi = initData['tths'][s]
        frm = initData['frms'][s]
        initSpot(k, s, t, etaRoi, tthRoi, frm, params, outPath, fname1, fname2)
    
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
            fname1 = os.path.join(dataDir,fnames[0])
            fname2 = os.path.join(dataDir,fnames[1])

            print('')
            print(f'Scan {t}, Spot:', end=" ")
            for k, s in enumerate(spotInds): 
                print(f'{k}', end=" ")
                # Save input file
                inputFile = os.path.join(outPath,'inputs',f'inputs_{k}.pkl')
                with open(inputFile, 'wb') as f:
                    pickle.dump((k,s,t,params,outPath,fname1,fname2), f)
                # Write job .sh file
                job_script = job_script_template.format(k=k,inputFile=inputFile)
                script_filename = os.path.join(outPath,'jobs',f'job_{k}.sh')
                with open(script_filename, 'w') as f:
                    f.write(job_script)
                # Submit job
                os.system(f"qsub {script_filename}")
 
            # Wait for my jobs to finish
            job_ids = collect_job_ids()
            wait_for_jobs(job_ids)
                
            i += 1
            t += 1

    
def compAvgParams(track):   
    J = len(track) 
    avgFWHMeta = 0
    avgFWHMtth = 0
    avgEta = 0
    avgTth = 0
    for j in range(J):
        avgFWHMeta += track[j]['p'][3]/J
        avgFWHMtth += track[j]['p'][4]/J
        avgEta += track[j]['eta']/J
        avgTth += track[j]['tth']/J
    
    return avgFWHMeta,avgFWHMtth,avgEta,avgTth
        