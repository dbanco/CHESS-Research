#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 10:14:50 2024

@author: dbanco
"""

import numpy as np
import scipy as sp
import yaml
import h5py
import hdf5plugin
from hexrd import imageseries

def read_yaml(file_path):
    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)
    return data

def loadYamlData(params,tth,eta):
    yamlFile = params['yamlFile']
    if params['detector'] == 'dexela':
        return loadYamlDataDexela(yamlFile,tth,eta)
    elif params['detector'] == 'eiger':
        return loadYamlDataEiger(yamlFile)
    elif params['detector'] == 'eiger_sim':
        return loadYamlDataEiger(yamlFile)

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

def loadYamlDataEiger(yamlFile):
    yamlData = read_yaml(yamlFile)
    ff_trans = yamlData['detectors']['eiger']['transform']['translation']
    detectDist = -ff_trans[2]
    mmPerPixel = yamlData['detectors']['eiger']['pixels']['size'][0]

    return detectDist, mmPerPixel, ff_trans

def load_eiger(fnames,params,t,detectSize=(4362,4148)):
    ims = imageseries.open(fnames[0], format='eiger-stream-v1')

    img = ims[t,:,:].copy()
    img[img>4294000000] = 0

    imSize = params['imSize']
    yamlFile = params['yamlFile']

    # Pad image
    bpad = np.zeros(imSize)
    center = (imSize[0]/2, imSize[1]/2)
    
    # Shift each panel
    detectDist, mmPerPixel, ff_trans = loadYamlDataEiger(yamlFile)
    ff1_tx = ff_trans[0]
    ff1_ty = ff_trans[1]
    ff1_xshift = int(round(ff1_tx/mmPerPixel))
    ff1_yshift = int(round(ff1_ty/mmPerPixel))
    
    # negative sign on y shift because rows increase downwards
    ff1r1 = int(center[0]-detectSize[0]/2-ff1_yshift)
    ff1r2 = int(center[0]+detectSize[0]/2-ff1_yshift)
    
    ff1c1 = int(center[1]-detectSize[1]/2+ff1_xshift)
    ff1c2 = int(center[1]+detectSize[1]/2+ff1_xshift)
        
    bpad[ff1r1:ff1r2,ff1c1:ff1c2] = img

    return bpad

def load_eiger_sim(fnames,params,t,detectSize=(4362,4148)):
    simData = np.load(fnames[0])       
    shp = simData['shape']
    img = np.zeros((shp[0],shp[1]))
    
    rowD = simData[f'{t}_row']
    colD = simData[f'{t}_col']
    datD = simData[f'{t}_data']
    
    for i in range(len(rowD)):
        img[rowD[i],colD[i]] = datD[i]

    imSize = params['imSize']
    yamlFile = params['yamlFile']

    # Pad image
    bpad = np.zeros(imSize)
    center = (imSize[0]/2, imSize[1]/2)
    
    # Shift each panel
    detectDist, mmPerPixel, ff_trans = loadYamlDataEiger(yamlFile)
    ff1_tx = ff_trans[0]
    ff1_ty = ff_trans[1]
    ff1_xshift = int(round(ff1_tx/mmPerPixel))
    ff1_yshift = int(round(ff1_ty/mmPerPixel))
    
    # negative sign on y shift because rows increase downwards
    ff1r1 = int(center[0]-detectSize[0]/2-ff1_yshift)
    ff1r2 = int(center[0]+detectSize[0]/2-ff1_yshift)
    
    ff1c1 = int(center[1]-detectSize[1]/2+ff1_xshift)
    ff1c2 = int(center[1]+detectSize[1]/2+ff1_xshift)
        
    bpad[ff1r1:ff1r2,ff1c1:ff1c2] = img

    return bpad    

def load_eiger_sim(fnames,params,t,detectSize=(4362,4148)):
    simData = np.load(fnames[0])       
    shp = simData['shape']
    img = np.zeros((shp[0],shp[1]))
    
    rowD = simData[f'{t}_row']
    colD = simData[f'{t}_col']
    datD = simData[f'{t}_data']
    
    for i in range(len(rowD)):
        img[rowD[i],colD[i]] = datD[i]

    imSize = params['imSize']
    yamlFile = params['yamlFile']

    # Pad image
    bpad = np.zeros(imSize)
    center = (imSize[0]/2, imSize[1]/2)
    
    # Shift each panel
    detectDist, mmPerPixel, ff_trans = loadYamlDataEiger(yamlFile)
    ff1_tx = ff_trans[0]
    ff1_ty = ff_trans[1]
    ff1_xshift = int(round(ff1_tx/mmPerPixel))
    ff1_yshift = int(round(ff1_ty/mmPerPixel))
    
    # negative sign on y shift because rows increase downwards
    ff1r1 = int(center[0]-detectSize[0]/2-ff1_yshift)
    ff1r2 = int(center[0]+detectSize[0]/2-ff1_yshift)
    
    ff1c1 = int(center[1]-detectSize[1]/2+ff1_xshift)
    ff1c2 = int(center[1]+detectSize[1]/2+ff1_xshift)
        
    bpad[ff1r1:ff1r2,ff1c1:ff1c2] = img

    return bpad  
    
def load_dex(fnames,params,t,t2=None,dexSize=(3888,3072)):
    if t2 is None:
        with h5py.File(fnames[0], 'r') as file1, h5py.File(fnames[1], 'r') as file2:
            img1 = file1['/imageseries/images'][t, :, :]
            img2 = file2['/imageseries/images'][t, :, :]
    else:
        with h5py.File(fnames[0], 'r') as file1, h5py.File(fnames[1], 'r') as file2:
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
    ff1_pix, ff2_pix = panelPixelsDex(ff_trans,mmPerPixel,imSize)
    
    new_center = np.array([center[0] - y_cart[0], center[1] - x_cart[0]])
    roiShape = getROIshapeDex(x_cart, y_cart, ff1_pix, ff2_pix, center)
    
    Ainterp = bilinearInterpMatrix(roiShape,rad_dom,eta_dom,new_center)
    
    return Ainterp, new_center, x_cart, y_cart

def getInterpParamsEiger(tth,eta,params):
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    imSize = params['imSize']
    
    center = (imSize[0]/2,imSize[1]/2)
    detectDist, mmPerPixel, ff_trans = loadYamlDataEiger(yamlFile)
    
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    x_cart, y_cart = fetchCartesian(rad_dom,eta_dom,center)
    ff_pix = panelPixelsEiger(ff_trans,mmPerPixel,imSize)
    
    new_center = np.array([center[0] - y_cart[0], center[1] - x_cart[0]])
    roiShape = getROIshapeEiger(x_cart, y_cart, ff_pix)
    
    Ainterp = bilinearInterpMatrix(roiShape,rad_dom,eta_dom,new_center)
    
    return Ainterp, new_center, x_cart, y_cart


def loadPolarROI(fnames,tth,eta,frame,params):
    if params['detector'] == 'dexela':
        roi = loadDexPolarRoi(fnames,tth,eta,frame,params)
    elif params['detector'] == 'eiger':
        roi = loadEigerPolarRoi(fnames[0],tth,eta,frame,params)
    elif params['detector'] == 'eiger_sim':
        roi = loadEigerSimPolarRoi(fnames[0],tth,eta,frame,params)
    return roi
    
def loadDexPolarRoi(fnames,tth,eta,frame,params):
    # 0. Load params, YAML data
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans = loadYamlDataDexela(yamlFile,tth,eta)
    
    # 1. Construct rad, eta domain
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    
    # 2. Construct interpolation matrix
    Ainterp,new_center,x_cart,y_cart = getInterpParamsDexela(tth,eta,params)
    
    # 3. Load needed Cartesian ROI pixels
    ff1_pix, ff2_pix = panelPixelsDex(ff_trans,mmPerPixel)
    roi = loadDexPanelROI(x_cart,y_cart,ff1_pix,ff2_pix,fnames,frame,params)
    
    # 4. Apply interpolation matrix to Cartesian pixels get Polar values
    roi_polar_vec = Ainterp.dot(roi.flatten())
    
    # 5. Reshape and output roi
    roi_polar = np.reshape(roi_polar_vec,(len(rad_dom),len(eta_dom)))
    
    return roi_polar

def loadEigerPolarRoi(fname,tth,eta,frame,params):
    # 0. Load params, YAML data
    roiSize = params['roiSize']
    imSize = params['imSize']
    detectDist, mmPerPixel, ff_trans = loadYamlData(params,tth,eta)
    
    # 1. Construct rad, eta domain
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    
    # 2. Construct interpolation matrix
    Ainterp,new_center,x_cart,y_cart = getInterpParamsEiger(tth,eta,params)
    
    # 3. Load needed Cartesian ROI pixels
    ff1_pix = panelPixelsEiger(ff_trans,mmPerPixel,imSize)
    roi = loadEigerPanelROI(x_cart,y_cart,ff1_pix,fname,frame)
    
    # 4. Apply interpolation matrix to Cartesian pixels get Polar values
    roi_polar_vec = Ainterp.dot(roi.flatten())
    
    # 5. Reshape and output roi
    roi_polar = np.reshape(roi_polar_vec,(len(rad_dom),len(eta_dom)))
    
    return roi_polar

def loadEigerSimPolarRoi(fname,tth,eta,frame,params):
    # 0. Load params, YAML data
    roiSize = params['roiSize']
    imSize = params['imSize']
    detectDist, mmPerPixel, ff_trans = loadYamlData(params,tth,eta)
    
    # 1. Construct rad, eta domain
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    
    # 2. Construct interpolation matrix
    Ainterp,new_center,x_cart,y_cart = getInterpParamsEiger(tth,eta,params)
    
    # 3. Load needed Cartesian ROI pixels
    ff1_pix = panelPixelsEiger(ff_trans,mmPerPixel,imSize)
    roi = loadEigerSimPanelROI(x_cart,y_cart,ff1_pix,fname,frame)
    
    # 4. Apply interpolation matrix to Cartesian pixels get Polar values
    roi_polar_vec = Ainterp.dot(roi.flatten())
    
    # 5. Reshape and output roi
    roi_polar = np.reshape(roi_polar_vec,(len(rad_dom),len(eta_dom)))
    
    return roi_polar

def loadEigerPolarRoiArray(fname,tth,eta,frames,params):
    # 0. Load params, YAML data
    imSize = params['imSize']
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans = loadYamlData(params,tth,eta)
    
    # 1. Construct rad, eta domain
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    
    # 2. Construct interpolation matrix
    Ainterp,new_center,x_cart,y_cart = getInterpParamsEiger(tth,eta,params)
    
    # 3. Load needed Cartesian ROI pixels
    ff1_pix = panelPixelsEiger(ff_trans,mmPerPixel,imSize)
    roiArray = loadEigerPanelROIArray(x_cart,y_cart,ff1_pix,fname,frames)
    roi_polar_array = np.zeros((frames[1]-frames[0],len(rad_dom),len(eta_dom)))
    # 4. Apply interpolation matrix to Cartesian pixels get Polar values
    for i in range(roiArray.shape[0]):
        print(i)
        roi_polar_vec = Ainterp.dot(roiArray[i,:,:].flatten())
        roi_polar = np.reshape(roi_polar_vec,(len(rad_dom),len(eta_dom)))
        roi_polar_array[i,:,:] = roi_polar
    
    return roi_polar_array

def loadDexPolarRoiPrecomp(fnames,yamlFile,tth,eta,frame,\
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
    ff1_pix, ff2_pix = panelPixelsDex(ff_trans,mmPerPixel)
    roi = loadDexPanelROI(x_cart,y_cart,ff1_pix,ff2_pix,fnames,frame,imSize)
    
    # 4. Apply interpolation matrix to Cartesian pixels get Polar values
    roi_polar_vec = Ainterp.dot(roi.flatten())
    
    # 5. Reshape and output roi
    roi_polar = np.reshape(roi_polar_vec,(len(rad_dom),len(eta_dom)))
    
    return roi_polar

def loadDexPanelROI(x_cart,y_cart,ff1_pix,ff2_pix,fnames,frame,params,dexShape=(3888,3072)):
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
        with h5py.File(fnames[1], 'r') as file:
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
        with h5py.File(fnames[0], 'r') as file:
            img = file['/imageseries/images'][frame, y_pan[0]:y_pan[1],x_pan[0]:x_pan[1]]
            img = np.fliplr(img)
      
    return img

def loadEigerPanelROI(x_cart,y_cart,ff1_pix,fname,frame):

    if x_cart[0] < ff1_pix[0]: x_cart[0] = ff1_pix[0]
    if x_cart[1] > ff1_pix[1]: x_cart[1] = ff1_pix[1]
    if y_cart[0] < ff1_pix[2]: y_cart[0] = ff1_pix[2]
    if y_cart[1] > ff1_pix[3]: y_cart[1] = ff1_pix[3]
    x_pan = x_cart - ff1_pix[0]
    y_pan = y_cart - ff1_pix[2]
    
    ims = imageseries.open(fname, format='eiger-stream-v1')
    img = ims[frame, y_pan[0]:y_pan[1],x_pan[0]:x_pan[1]].copy()
    img[img>4294000000] = 0
    return img

def loadEigerSimPanelROI(x_cart,y_cart,ff1_pix,fname,frame):

    if x_cart[0] < ff1_pix[0]: x_cart[0] = ff1_pix[0]
    if x_cart[1] > ff1_pix[1]: x_cart[1] = ff1_pix[1]
    if y_cart[0] < ff1_pix[2]: y_cart[0] = ff1_pix[2]
    if y_cart[1] > ff1_pix[3]: y_cart[1] = ff1_pix[3]
    x_pan = x_cart - ff1_pix[0]
    y_pan = y_cart - ff1_pix[2]
    
    simData = np.load(fname)       
    shp = simData['shape']
    imgFull = np.zeros((shp[0],shp[1]))

    frame = frame - 2
    rowD = simData[f'{frame}_row']
    colD = simData[f'{frame}_col']
    datD = simData[f'{frame}_data']

    for i in range(len(rowD)):
        imgFull[rowD[i],colD[i]] = datD[i]
    
    img = imgFull[y_pan[0]:y_pan[1],x_pan[0]:x_pan[1]]
    
    return img

def loadEigerPanelROIArray(x_cart,y_cart,ff1_pix,fname,frames):

    if x_cart[0] < ff1_pix[0]: x_cart[0] = ff1_pix[0]
    if x_cart[1] > ff1_pix[1]: x_cart[1] = ff1_pix[1]
    if y_cart[0] < ff1_pix[2]: y_cart[0] = ff1_pix[2]
    if y_cart[1] > ff1_pix[3]: y_cart[1] = ff1_pix[3]
    x_pan = x_cart - ff1_pix[0]
    y_pan = y_cart - ff1_pix[2]
    
    ims = imageseries.open(fname, format='eiger-stream-v1')
    imgArray = np.zeros((frames[1]-frames[0],y_pan[1]-y_pan[0],x_pan[1]-x_pan[0]))
    for i in np.arange(frames[0],frames[1]):
        print(i)
        imgArray[i-frames[0],:,:] = ims[i, y_pan[0]:y_pan[1],x_pan[0]:x_pan[1]].copy()
    imgArray[imgArray>4294000000] = 0
    
    return imgArray

def getROIshapeDex(x_cart,y_cart,ff1_pix,ff2_pix,center,dexShape=(3888,3072)):
    
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

def getROIshapeEiger(x_cart,y_cart,ff1_pix,eigShape=(4362,4148)):
    
    if x_cart[0] < ff1_pix[0]: x_cart[0] = ff1_pix[0]
    if x_cart[1] > ff1_pix[1]: x_cart[1] = ff1_pix[1]
    if y_cart[0] < ff1_pix[2]: y_cart[0] = ff1_pix[2]
    if y_cart[1] > ff1_pix[3]: y_cart[1] = ff1_pix[3]
    x_pan = x_cart - ff1_pix[0]
    y_pan = y_cart - ff1_pix[2]

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

def panelPixelsDex(ff_trans,mmPerPixel,imSize=(4888,7300),dexShape=(3888,3072)):
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

def panelPixelsEiger(ff_trans,mmPerPixel,imSize=(5000,5000),detectShape=(4362,4148)):
    center = (imSize[0]/2,imSize[1]/2)
    ff1_tx = ff_trans[0]
    ff1_ty = ff_trans[1]
    
    ff1_xshift = int(round(ff1_tx/mmPerPixel))
    ff1_yshift = int(round(ff1_ty/mmPerPixel))
    
    # negative sign on y shift because rows increase downwards
    ff1r1 = int(center[0]-detectShape[0]/2-ff1_yshift)
    ff1r2 = int(center[0]+detectShape[0]/2-ff1_yshift)
    
    ff1c1 = int(center[1]-detectShape[1]/2+ff1_xshift)
    ff1c2 = int(center[1]+detectShape[1]/2+ff1_xshift)
    
    ff1_pixels = [ff1c1,ff1c2,ff1r1,ff1r2]
    return ff1_pixels

def polarDomain(detectDist,mmPerPixel,tth,eta,roi_size):
    r = np.round(detectDist*np.tan(tth)/mmPerPixel)
    r1 = r - roi_size[0]//2
    r2 = r + roi_size[0]//2

    rad_domain = np.arange(r1,r2,1)
    
    deta = 1/r2
    eta1 = eta - roi_size[1]/2*deta
    eta2 = eta + roi_size[1]/2*deta
    eta_domain = np.linspace(eta1,eta2,roi_size[1])
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
