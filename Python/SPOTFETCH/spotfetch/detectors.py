#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spotfetch.detectors

Description:
Models reading data from x-ray detectors

Supported detectors:
    dexela
    eiger
    eiger_sim
    
Created on: Fri Nov 15 23:00:28 2024
Author: Daniel Banco
email: dpqb10@gmail.com
"""
import numpy as np
import scipy as sp
import yaml
import h5py
import os
import glob
from hexrd import imageseries
from hexrd import transforms
import spotfetch as sf

DEX_SHAPE = (3888, 3072)
EIG_SHAPE = (4362, 4148)

def read_yaml(file_path):
    """
    Reads a YAML file and returns the parsed data.

    Parameters:
    -----------
    file_path : str
        Path to the YAML file.

    Returns:
    --------
    dict
        Parsed YAML data.
    """
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

def loadYamlData(params, tth=None, eta=None):
    """
    Loads YAML data based on detector type and parameters.

    Parameters:
    -----------
    params : dict
        Parameters including detector type and YAML file path.
    tth : float, optional
        Theta angle for the detector (used for Dexela).
    eta : float, optional
        Eta angle for the detector (used for Dexela).

    Returns:
    --------
    tuple
        Detector distance, pixel size, translation, and tilt.
    """
    yamlFile = params['yamlFile']
    detector = params['detector']

    if detector == 'dexela':
        if tth is None:
            raise ValueError("Theta (tth) must be provided for Dexela detector.")
        return loadYamlDataDexela(yamlFile, tth, eta)
    elif detector in ['eiger', 'eiger_sim']:
        return loadYamlDataEiger(yamlFile)
    else:
        raise ValueError(f"Unsupported detector type: {detector}")
        
def loadYamlDataDexela(yamlFile, tth, eta):
    """
    Loads YAML data for the Dexela detector.

    Parameters:
    -----------
    yamlFile : str
        Path to the YAML file.
    tth : float
        Theta angle for the detector.
    eta : float
        Eta angle for the detector.

    Returns:
    --------
    tuple
        Detector distance, pixel size, translation, and tilt.
    """
    yamlData = read_yaml(yamlFile)
    ff1_trans = yamlData['detectors']['ff1']['transform']['translation']
    ff2_trans = yamlData['detectors']['ff2']['transform']['translation']
    ff1_tilt = yamlData['detectors']['ff1']['transform']['tilt']
    ff2_tilt = yamlData['detectors']['ff2']['transform']['tilt']

    ff_trans = [ff1_trans[0], ff1_trans[1], ff2_trans[0], ff2_trans[1]]
    ff_tilt = [ff1_tilt, ff2_tilt]

    if abs(eta) > (np.pi / 2):
        detectDist = -ff2_trans[2]
        mmPerPixel = yamlData['detectors']['ff2']['pixels']['size'][0]
    else:
        detectDist = -ff1_trans[2]
        mmPerPixel = yamlData['detectors']['ff1']['pixels']['size'][0]

    return detectDist, mmPerPixel, ff_trans, ff_tilt

def loadYamlDataEiger(yamlFile):
    """
    Loads YAML data for the Eiger detector.

    Parameters:
    -----------
    yamlFile : str
        Path to the YAML file.

    Returns:
    --------
    tuple
        Detector distance, pixel size, translation, and tilt.
    """
    yamlData = read_yaml(yamlFile)
    trans = yamlData['detectors']['eiger']['transform']['translation']
    tilt = yamlData['detectors']['eiger']['transform']['tilt']
    detectDist = -trans[2]
    mmPerPixel = yamlData['detectors']['eiger']['pixels']['size'][0]

    return detectDist, mmPerPixel, trans, tilt

def loadImg(fnames, params, frame):
    """
    Loads an image based on the detector type.

    Parameters:
    -----------
    fnames : list
        File names for image files (one file per detector panel)
    params : dict
        Parameters including detector type.
    frame : int
        Frame index to load.

    Returns:
    --------
    ndarray
        Loaded image.
    """
    detector = params['detector']

    if detector == 'dexela':
        return loadDexImg(fnames, params, frame)
    elif detector == 'eiger':
        return loadEigerImg(fnames, params, frame)
    elif detector == 'eiger_sim':
        return loadEigerSimImg(fnames, params, frame)
    else:
        raise ValueError(f"Unsupported detector type: {detector}")
        
def loadDexImg(fnames, params, frame_i, dexSize=DEX_SHAPE):
    """
    Loads an image for the Dexela detector and aligns its panels.

    Parameters:
    -----------
    fnames : list
        List of file names (one for each panel).
    params : dict
        Parameters including image size and YAML file path.
    frame_i : int
        Frame index to load.
    dexSize : tuple, optional
        Size of each Dexela panel, default is DEX_SHAPE.

    Returns:
    --------
    ndarray
        Padded and aligned image.
    """
    with h5py.File(fnames[0], 'r') as file1, h5py.File(fnames[1], 'r') as file2:
        img1 = file1['/imageseries/images'][frame_i, :, :]
        img2 = file2['/imageseries/images'][frame_i, :, :]

    imSize = params['imSize']
    yamlFile = params['yamlFile']
    bpad = np.zeros(imSize)
    center = (imSize[0] / 2, imSize[1] / 2)

    yamlData = read_yaml(yamlFile)
    mmPerPixel = yamlData['detectors']['ff2']['pixels']['size'][0]
    ff1_trans = yamlData['detectors']['ff1']['transform']['translation']
    ff2_trans = yamlData['detectors']['ff2']['transform']['translation']

    shifts = {
        'ff1_x': int(round(ff1_trans[0] / mmPerPixel)),
        'ff1_y': int(round(ff1_trans[1] / mmPerPixel)),
        'ff2_x': int(round(ff2_trans[0] / mmPerPixel)),
        'ff2_y': int(round(ff2_trans[1] / mmPerPixel)),
    }

    # Define panel coordinates
    coords = {
        'ff1': {
            'r1': int(center[0] - dexSize[0] / 2 - shifts['ff1_y']),
            'r2': int(center[0] + dexSize[0] / 2 - shifts['ff1_y']),
            'c1': int(center[1] - dexSize[1] / 2 + shifts['ff1_x']),
            'c2': int(center[1] + dexSize[1] / 2 + shifts['ff1_x']),
        },
        'ff2': {
            'r1': int(center[0] - dexSize[0] / 2 - shifts['ff2_y']),
            'r2': int(center[0] + dexSize[0] / 2 - shifts['ff2_y']),
            'c1': int(center[1] - dexSize[1] / 2 + shifts['ff2_x']),
            'c2': int(center[1] + dexSize[1] / 2 + shifts['ff2_x']),
        }
    }

    # Assign flipped images to padded array
    bpad[coords['ff2']['r1']:coords['ff2']['r2'], coords['ff2']['c1']:coords['ff2']['c2']] = np.flipud(img2)
    bpad[coords['ff1']['r1']:coords['ff1']['r2'], coords['ff1']['c1']:coords['ff1']['c2']] = np.fliplr(img1)

    return bpad

def loadEigerImg(fnames, params, frame, detectSize=EIG_SHAPE):
    """
    Loads and aligns an image for the Eiger detector.

    Parameters:
    -----------
    fnames : str
        File name of the Eiger image series.
    params : dict
        Parameters including image size and YAML file path.
    frame : int
        Frame index to load.
    detectSize : tuple, optional
        Size of the Eiger detector, default is EIG_SHAPE.

    Returns:
    --------
    ndarray
        Padded and aligned image.
    """
    # 0. Load params, YAML data
    imSize = params['imSize']
    yamlFile = params['yamlFile']
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlDataEiger(yamlFile)

    ims = imageseries.open(fnames, format='eiger-stream-v1')
    img = ims[frame, :, :].copy()
    img[img > 4294000000] = 0

    # Pad image
    bpad = np.zeros(imSize)
    center = (imSize[0] / 2, imSize[1] / 2)

    # Shift each panel
    ff1_tx = ff_trans[0]
    ff1_ty = ff_trans[1]
    ff1_xshift = int(round(ff1_tx / mmPerPixel))
    ff1_yshift = int(round(ff1_ty / mmPerPixel))

    # Negative sign on y shift because rows increase downwards
    ff1r1 = int(center[0] - detectSize[0] / 2 - ff1_yshift)
    ff1r2 = int(center[0] + detectSize[0] / 2 - ff1_yshift)

    ff1c1 = int(center[1] - detectSize[1] / 2 + ff1_xshift)
    ff1c2 = int(center[1] + detectSize[1] / 2 + ff1_xshift)

    bpad[ff1r1:ff1r2, ff1c1:ff1c2] = img

    return bpad

def loadEigerSimImg(fnames, params, frame_i, detectSize=EIG_SHAPE):
    """
    Loads and aligns a simulated image for the Eiger detector.

    Parameters:
    -----------
    fnames : str
        File name of the simulated data (.npz format).
    params : dict
        Parameters including image size and YAML file path.
    frame_i : int
        Frame index to load.
    detectSize : tuple, optional
        Size of the Eiger detector, default is EIG_SHAPE.

    Returns:
    --------
    ndarray
        Padded and aligned image.
    """
    simData = np.load(fnames)
    shp = simData['shape']
    img = np.zeros((shp[0], shp[1]))

    rowD = simData[f'{frame_i}_row']
    colD = simData[f'{frame_i}_col']
    datD = simData[f'{frame_i}_data']

    for i in range(len(rowD)):
        img[rowD[i], colD[i]] = datD[i]

    imSize = params['imSize']
    yamlFile = params['yamlFile']

    # Pad image
    bpad = np.zeros(imSize)
    center = (imSize[0] / 2, imSize[1] / 2)

    # Shift each panel
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlDataEiger(yamlFile)
    ff1_tx = ff_trans[0]
    ff1_ty = ff_trans[1]
    ff1_xshift = int(round(ff1_tx / mmPerPixel))
    ff1_yshift = int(round(ff1_ty / mmPerPixel))

    # Negative sign on y shift because rows increase downwards
    ff1r1 = int(center[0] - detectSize[0] / 2 - ff1_yshift)
    ff1r2 = int(center[0] + detectSize[0] / 2 - ff1_yshift)

    ff1c1 = int(center[1] - detectSize[1] / 2 + ff1_xshift)
    ff1c2 = int(center[1] + detectSize[1] / 2 + ff1_xshift)

    bpad[ff1r1:ff1r2, ff1c1:ff1c2] = img

    return bpad

def getInterpParams(tth, eta, params):
    """
    Retrieves interpolation parameters for the detector.

    Parameters:
    -----------
    tth : float
        Two-theta value.
    eta : float
        Eta value.
    params : dict
        Parameters including detector type and YAML file path.

    Returns:
    --------
    tuple
        Interpolation matrix, new center, x_cart, and y_cart coordinates.
    """
    if params['detector'] == 'dexela':
        return getInterpParamsDexela(tth, eta, params)
    elif params['detector'] == 'eiger':
        return getInterpParamsEiger(tth, eta, params)
    elif params['detector'] == 'eiger_sim':
        return getInterpParamsEiger(tth, eta, params)

def getInterpParamsDexela(tth, eta, params):
    """
    Retrieves interpolation parameters for the Dexela detector.

    Parameters:
    -----------
    tth : float
        Two-theta value.
    eta : float
        Eta value.
    params : dict
        Parameters including detector type and YAML file path.

    Returns:
    --------
    tuple
        Interpolation matrix, new center, x_cart, and y_cart coordinates.
    """
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    imSize = params['imSize']

    center = (imSize[0] / 2, imSize[1] / 2)
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlDataDexela(yamlFile, tth, eta)

    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    x_cart, y_cart = fetchCartesian(rad_dom, eta_dom, center)
    ff1_pix, ff2_pix = panelPixelsDex(ff_trans, mmPerPixel, imSize)

    new_center = np.array([center[0] - y_cart[0], center[1] - x_cart[0]])
    roiShape = getROIshapeDex(x_cart, y_cart, ff1_pix, ff2_pix, center)

    Ainterp = bilinearInterpMatrix(roiShape, rad_dom, eta_dom, new_center, detectDist)

    return Ainterp, new_center, x_cart, y_cart

def getInterpParamsEiger(tth, eta, params):
    """
    Retrieves interpolation parameters for the Eiger detector.

    Parameters:
    -----------
    tth : float
        Two-theta value.
    eta : float
        Eta value.
    params : dict
        Parameters including YAML file path, ROI size, and image size.

    Returns:
    --------
    tuple
        Interpolation matrix, new center, x_cart, and y_cart coordinates.
    """
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    imSize = params['imSize']

    center = ((imSize[0] - 1) / 2, (imSize[1] - 1) / 2)
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlDataEiger(yamlFile)

    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)

    # Get interpolation matrix
    x_cart, y_cart = fetchCartesian(rad_dom, eta_dom, center)
    ff_pix = panelPixelsEiger(ff_trans, mmPerPixel, imSize)

    new_center = np.array([center[0] - y_cart[0], center[1] - x_cart[0]])
    x_pan, y_pan = getEigerPixels(x_cart, y_cart, ff_pix)
    roiShape = [y_pan[1]-y_pan[0],x_pan[1]-x_pan[0]]

    Ainterp = bilinearInterpMatrix(roiShape, rad_dom, eta_dom, new_center, detectDist)

    return Ainterp, new_center, x_cart, y_cart

def loadROI(dataPath, scan, frame, etaRoi, tthRoi, params):
    """
    Loads the region of interest (ROI) for the given detector.

    Parameters:
    -----------
    dataPath : str
        Path to the data.
    scan : int
        Scan number.
    frame : int
        Frame index.
    etaRoi : float
        Eta range of interest.
    tthRoi : float
        Two-theta range of interest.
    params : dict
        Detector parameters.

    Returns:
    --------
    ndarray
        Loaded polar ROI.
    """
    if os.path.isdir(dataPath):
        fnames = sf.timeToFile(scan, dataPath)
        isFile = False
    else:
        fnames = dataPath
        isFile = True
        if params['detector'] == 'eiger_sim':
            fnames = fnames.replace('*', '{scan}').format(scan=scan)
            roi = loadPolarROI(fnames, tthRoi, etaRoi, frame, params)
            return roi

    # Load ROI
    if isFile:
        template = dataPath.format(num2=scan)
        fnames = glob.glob(template)
        roi = loadPolarROI(fnames[0], tthRoi, etaRoi, frame, params)
    else:
        roi = loadPolarROI(fnames, tthRoi, etaRoi, frame, params)

    return roi

def loadPolarROI(fnames, tth, eta, frame, params):
    """
    Loads the polar ROI for the specified detector type.

    Parameters:
    -----------
    fnames : str or list
        File name(s) of the image data.
    tth : float
        Two-theta value.
    eta : float
        Eta value.
    frame : int
        Frame index.
    params : dict
        Detector parameters.

    Returns:
    --------
    ndarray
        Loaded polar ROI.
    """
    if params['detector'] == 'dexela':
        roi = loadDexPolarRoi(fnames, tth, eta, frame, params)
    elif params['detector'] == 'eiger':
        roi = loadEigerPolarRoi(fnames, tth, eta, frame, params)
    elif params['detector'] == 'eiger_sim':
        roi = loadEigerPolarRoi(fnames, tth, eta, frame, params)
    return roi

def loadPolarROI3D(fnames, tth, eta, frame, params):
    """
    Loads a 3D polar ROI for the specified detector type.

    Parameters:
    -----------
    fnames : str or list
        File name(s) of the image data.
    tth : float
        Two-theta value.
    eta : float
        Eta value.
    frame : int
        Central frame index.
    params : dict
        Detector parameters.

    Returns:
    --------
    ndarray
        Loaded 3D polar ROI.
    """
    if params['detector'] == 'dexela':
        roi3D = loadDexPolarRoi3D(fnames, tth, eta, frame, params)
    elif params['detector'] == 'eiger':
        roi3D = loadEigerPolarRoi3D(fnames, tth, eta, frame, params)
    elif params['detector'] == 'eiger_sim':
        roi3D = loadEigerSimPolarRoi3D(fnames, tth, eta, frame, params)
    return roi3D

def loadDexPolarRoi(fnames, tth, eta, frame, params):
    """
    Loads the polar ROI for the Dexela detector.

    Parameters:
    -----------
    fnames : str or list
        File name(s) of the image data.
    tth : float
        Two-theta value.
    eta : float
        Eta value.
    frame : int
        Frame index.
    params : dict
        Detector parameters.

    Returns:
    --------
    ndarray
        Loaded polar ROI.
    """
    # Load parameters and YAML data
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlDataDexela(yamlFile, tth, eta)

    # Construct radial and eta domains
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)

    # Construct interpolation matrix
    Ainterp, new_center, x_cart, y_cart = getInterpParamsDexela(tth, eta, params)

    # Load Cartesian ROI pixels
    ff1_pix, ff2_pix = panelPixelsDex(ff_trans, mmPerPixel)
    roi = loadDexPanelROI(x_cart, y_cart, ff1_pix, ff2_pix, fnames, frame, params)

    # Apply interpolation matrix to Cartesian pixels to get polar values
    roi_polar_vec = Ainterp.dot(roi.flatten())

    # Reshape and return ROI
    roi_polar = np.reshape(roi_polar_vec, (len(rad_dom), len(eta_dom)))

    return roi_polar

def loadDexPolarRoi3D(fnames, tth, eta, frame, params):
    """
    Loads a 3D polar ROI for the Dexela detector.

    Parameters:
    -----------
    fnames : str or list
        File name(s) of the image data.
    tth : float
        Two-theta value.
    eta : float
        Eta value.
    frame : int
        Central frame index.
    params : dict
        Detector parameters.

    Returns:
    --------
    ndarray
        Loaded 3D polar ROI.
    """
    # Load parameters and YAML data
    yamlFile = params['yamlFile']
    roiSize = params['roiSize']
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlDataDexela(yamlFile, tth, eta)
    dome = (roiSize[2] - 1) / 2
    frmRange = sf.wrapFrame(np.arange(frame - dome, frame + dome + 1))

    # Construct radial and eta domains
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)

    # Construct interpolation matrix
    Ainterp, new_center, x_cart, y_cart = getInterpParamsDexela(tth, eta, params)

    # Load Cartesian ROI pixels for each frame
    ff1_pix, ff2_pix = panelPixelsDex(ff_trans, mmPerPixel)
    roi3D = np.zeros(roiSize)
    for i, frm in enumerate(frmRange):
        roi = loadDexPanelROI(x_cart, y_cart, ff1_pix, ff2_pix, fnames, int(frm), params)
        roi_polar_vec = Ainterp.dot(roi.flatten())
        roi_polar = np.reshape(roi_polar_vec, roiSize[:2])
        roi3D[:, :, i] = roi_polar

    return roi3D

def loadEigerPolarRoi(fname, tth, eta, frame, params):
    """
    Loads a polar region of interest (ROI) for the Eiger detector.

    Parameters:
    -----------
    fname : str
        File name of the data.
    tth : float
        Two-theta value.
    eta : float
        Eta value.
    frame : int
        Frame index to load.
    params : dict
        Parameters including ROI size, image size, and YAML file path.

    Returns:
    --------
    ndarray
        Reshaped polar ROI.
    """
    # 0. Load params, YAML data
    roiSize = params['roiSize']
    imSize = params['imSize']
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlData(params)
    
    # 1. Construct rad, eta domain
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    
    # 2. Construct interpolation matrix
    Ainterp, new_center, x_cart, y_cart = getInterpParamsEiger(tth, eta, params)
    
    # 3. Load needed Cartesian ROI pixels
    ff1_pix = panelPixelsEiger(ff_trans, mmPerPixel, imSize)
    if params['detector'] == 'eiger':
        roi = loadEigerPanelROI(x_cart, y_cart, ff1_pix, fname, frame)
    elif params['detector'] == 'eiger_sim':
        roi = loadEigerSimPanelROI(x_cart, y_cart, ff1_pix, fname, frame)
    
    # 4. Apply interpolation matrix to Cartesian pixels get Polar values
    roi_polar_vec = Ainterp.dot(roi.flatten())
    
    # 5. Reshape and output roi
    roi_polar = np.reshape(roi_polar_vec, (len(rad_dom), len(eta_dom)))

    return roi_polar

def loadEigerPolarRoi3D(fname,tth,eta,frame,params):
    """
    Loads a 3D polar region of interest (ROI) for the Eiger detector.

    Parameters:
    -----------
    fname : str
        File name of the data.
    tth : float
        Two-theta value.
    eta : float
        Eta value.
    frame : int
        Frame index to load.
    params : dict
        Parameters including ROI size, image size, and YAML file path.

    Returns:
    --------
    ndarray
        3D polar ROI.
    """
    # 0. Load params, YAML data
    roiSize = params['roiSize']
    imSize = params['imSize']
    dome = (roiSize[2]-1)/2
    frmRange = sf.wrapFrame(np.arange(frame-dome,frame+dome+1))
    
    detectDist, mmPerPixel, ff_trans, ff_tilt = loadYamlData(params)
    
    # 1. Construct rad, eta domain
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tth, eta, roiSize)
    
    # 2. Construct interpolation matrix
    Ainterp,new_center,x_cart,y_cart = getInterpParamsEiger(tth,eta,params)
    
    ff1_pix = panelPixelsEiger(ff_trans,mmPerPixel,imSize)
    roi3D = np.zeros(roiSize)
    for i, frm in enumerate(frmRange):
        # 3. Load needed Cartesian ROI pixels
        if params['detector'] == 'eiger':
            roi = loadEigerPanelROI(x_cart, y_cart, ff1_pix, fname, frame)
        elif params['detector'] == 'eiger_sim':
            roi = loadEigerSimPanelROI(x_cart, y_cart, ff1_pix, fname, frame)
        # 4. Apply interpolation matrix to Cartesian pixels get Polar values
        roi_polar_vec = Ainterp.dot(roi.flatten())
        # 5. Reshape and output roi
        roi_polar = np.reshape(roi_polar_vec, roiSize[:2])
        roi3D[:, :, i] = roi_polar
    
    return roi3D

def loadDexPanelROI(x_cart, y_cart, ff1_pix, ff2_pix, fnames, frame, params, dexShape=DEX_SHAPE):
    """
    Loads a Region of Interest (ROI) for the Dexela detector in Cartesian coordinates.

    Parameters:
    -----------
    x_cart : list or ndarray
        X-coordinates of the ROI in the Cartesian system.
    y_cart : list or ndarray
        Y-coordinates of the ROI in the Cartesian system.
    ff1_pix : ndarray
        Flat-field parameters for the first panel (right panel).
    ff2_pix : ndarray
        Flat-field parameters for the second panel (left panel).
    fnames : list
        List of file names containing panel data.
    frame : int
        Frame index to load.
    params : dict
        Dictionary of detector parameters, including image size.
    dexShape : tuple, optional
        Shape of the Dexela detector, default is DEX_SHAPE.

    Returns:
    --------
    ndarray
        Loaded ROI image from the specified panel.
    """
    # Extract image size and calculate the center
    imSize = params['imSize']
    center = (imSize[0] / 2, imSize[1] / 2)

    # Ensure ROI boundaries do not exceed panel limits
    if x_cart[0] < ff2_pix[0]: x_cart[0] = ff2_pix[0]
    if x_cart[1] > ff1_pix[1]: x_cart[1] = ff1_pix[1]
    if y_cart[0] < ff2_pix[2]: y_cart[0] = ff2_pix[2]
    if y_cart[1] > ff2_pix[3]: y_cart[1] = ff2_pix[3]

    # Determine which panel the ROI belongs to
    if x_cart[0] < center[1]:  # Left panel
        # Adjust for panel offsets
        x_pan = x_cart - ff2_pix[0]
        y_pan = y_cart - ff2_pix[2]
        # Account for vertical flipping
        midLine = (dexShape[0] - 1) / 2
        flip0 = y_pan[0] + 2 * (midLine - y_pan[0])
        flip1 = y_pan[1] + 2 * (midLine - y_pan[1])
        y_pan[0] = min(flip0, flip1)
        y_pan[1] = max(flip0, flip1)
        # Load and flip the image vertically
        with h5py.File(fnames[1], 'r') as file:
            img = file['/imageseries/images'][frame, y_pan[0]:y_pan[1], x_pan[0]:x_pan[1]]
            img = np.flipud(img)
    elif x_cart[0] > center[1]:  # Right panel
        # Adjust for panel offsets
        x_pan = x_cart - ff1_pix[0]
        y_pan = y_cart - ff1_pix[2]
        # Account for horizontal flipping
        midLine = (dexShape[1] - 1) / 2
        flip0 = x_pan[0] + 2 * (midLine - x_pan[0])
        flip1 = x_pan[1] + 2 * (midLine - x_pan[1])
        x_pan[0] = min(flip0, flip1)
        x_pan[1] = max(flip0, flip1)
        # Load and flip the image horizontally
        with h5py.File(fnames[0], 'r') as file:
            img = file['/imageseries/images'][frame, y_pan[0]:y_pan[1], x_pan[0]:x_pan[1]]
            img = np.fliplr(img)

    return img

def loadEigerPanelROI(x_cart, y_cart, ff1_pix, fname, frame):
    """
    Loads a Region of Interest (ROI) for the Eiger detector in Cartesian coordinates.

    Parameters:
    -----------
    x_cart : list or ndarray
        X-coordinates of the ROI in the Cartesian system.
    y_cart : list or ndarray
        Y-coordinates of the ROI in the Cartesian system.
    ff1_pix : ndarray
        Flat-field parameters for the panel.
    fname : str
        File name of the image series.
    frame : int
        Frame index to load.

    Returns:
    --------
    ndarray
        Loaded ROI image.
    """
    # Compute panel-specific coordinates
    x_pan, y_pan = getEigerPixels(x_cart, y_cart, ff1_pix)

    # Load image data
    ims = imageseries.open(fname, format='eiger-stream-v1')
    img = ims[frame, y_pan[0]:y_pan[1], x_pan[0]:x_pan[1]].copy()

    # Mask invalid pixel values
    img[img > 4294000000] = 0
    return img

def loadEigerSimPanelROI(x_cart, y_cart, ff1_pix, fname, frame):
    """
    Loads a simulated Region of Interest (ROI) for the Eiger detector.

    Parameters:
    -----------
    x_cart : list or ndarray
        X-coordinates of the ROI in the Cartesian system.
    y_cart : list or ndarray
        Y-coordinates of the ROI in the Cartesian system.
    ff1_pix : ndarray
        Flat-field parameters for the panel.
    fname : str
        File name of the simulation data (NumPy file).
    frame : int
        Frame index to load.

    Returns:
    --------
    ndarray
        Loaded simulated ROI image.
    """
    # Compute panel-specific coordinates
    x_pan, y_pan = getEigerPixels(x_cart, y_cart, ff1_pix)

    # Load simulation data
    simData = np.load(fname)
    shp = simData['shape']
    imgFull = np.zeros((shp[0], shp[1]))

    # Load data for the specific frame
    frame = frame - 2
    rowD = simData[f'{frame}_row']
    colD = simData[f'{frame}_col']
    datD = simData[f'{frame}_data']

    # Populate the full image array
    for i in range(len(rowD)):
        imgFull[rowD[i], colD[i]] = datD[i]

    # Extract and return the ROI
    img = imgFull[y_pan[0]:y_pan[1], x_pan[0]:x_pan[1]].copy()
    return img

def loadEigerPanelROIArray(x_cart, y_cart, ff1_pix, fname, frames):
    """
    Loads a series of Regions of Interest (ROIs) for the Eiger detector across multiple frames.

    Parameters:
    -----------
    x_cart : list or ndarray
        X-coordinates of the ROI in the Cartesian system.
    y_cart : list or ndarray
        Y-coordinates of the ROI in the Cartesian system.
    ff1_pix : ndarray
        Flat-field parameters for the panel.
    fname : str
        File name of the image series.
    frames : tuple
        Start and end frame indices (inclusive).

    Returns:
    --------
    ndarray
        3D array containing ROIs for the specified frames.
    """
    # Compute panel-specific coordinates
    x_pan, y_pan = getEigerPixels(x_cart, y_cart, ff1_pix)

    # Load image data
    ims = imageseries.open(fname, format='eiger-stream-v1')
    imgArray = np.zeros((frames[1] - frames[0], y_pan[1] - y_pan[0], x_pan[1] - x_pan[0]))

    # Extract ROI for each frame
    for i in np.arange(frames[0], frames[1]):
        imgArray[i - frames[0], :, :] = ims[i, y_pan[0]:y_pan[1], x_pan[0]:x_pan[1]].copy()

    # Mask invalid pixel values
    imgArray[imgArray > 4294000000] = 0
    return imgArray

def loadEigerSimPanelROIArray(x_cart, y_cart, ff1_pix, fname, frames):
    """
    Loads a series of simulated Regions of Interest (ROIs) for the Eiger detector across multiple frames.

    Parameters:
    -----------
    x_cart : list or ndarray
        X-coordinates of the ROI in the Cartesian system.
    y_cart : list or ndarray
        Y-coordinates of the ROI in the Cartesian system.
    ff1_pix : ndarray
        Flat-field parameters for the panel.
    fname : str
        File name of the simulation data (NumPy file).
    frames : tuple
        Start and end frame indices (inclusive).

    Returns:
    --------
    ndarray
        3D array containing simulated ROIs for the specified frames.
    """
    # Compute panel-specific coordinates
    x_pan, y_pan = getEigerPixels(x_cart, y_cart, ff1_pix)

    # Load simulation data
    simData = np.load(fname)
    shp = simData['shape']
    imgFull = np.zeros((frames[1] - frames[0], shp[0], shp[1]))

    # Extract ROI for each frame in the specified range
    for frm in np.arange(frames[0], frames[1]):
        frm = frm - 2  # Adjust frame index
        rowD = simData[f'{frm}_row']
        colD = simData[f'{frm}_col']
        datD = simData[f'{frm}_data']

        # Populate the full image array for the frame
        for i in range(len(rowD)):
            imgFull[2 + frm - frames[0], rowD[i], colD[i]] = datD[i]

    # Extract and return the ROI for all frames
    return imgFull[:, y_pan[0]:y_pan[1], x_pan[0]:x_pan[1]].copy()

def getROIshapeDex(x_cart, y_cart, ff1_pix, ff2_pix, center, dexShape=DEX_SHAPE):
    """
    Calculates the shape of the ROI for the Dexela detector based on Cartesian coordinates.

    Parameters:
    -----------
    x_cart : list or ndarray
        X-coordinates of the ROI in the Cartesian system.
    y_cart : list or ndarray
        Y-coordinates of the ROI in the Cartesian system.
    ff1_pix : list or ndarray
        Flat-field parameters for the primary panel.
    ff2_pix : list or ndarray
        Flat-field parameters for the secondary panel.
    center : tuple
        Image center coordinates.
    dexShape : tuple, optional
        Shape of the Dexela detector panel (default is DEX_SHAPE).

    Returns:
    --------
    tuple
        Shape of the ROI (rows, columns).
    """
    # Ensure x-coordinates are within bounds
    if x_cart[0] < ff2_pix[0]: 
        x_cart[0] = ff2_pix[0]
    if x_cart[1] > ff1_pix[1]: 
        x_cart[1] = ff1_pix[1]

    # Determine the panel and process accordingly
    if x_cart[0] < center[1]:  # Panel 2
        if y_cart[0] < ff2_pix[2]: 
            y_cart[0] = ff2_pix[2]
        if y_cart[1] > ff2_pix[3]: 
            y_cart[1] = ff2_pix[3]
        x_pan = x_cart - ff2_pix[0]
        y_pan = y_cart - ff2_pix[2]

        # Account for vertical flipping
        midLine = (dexShape[0] - 1) / 2
        flip0 = y_pan[0] + 2 * (midLine - y_pan[0])
        flip1 = y_pan[1] + 2 * (midLine - y_pan[1])
        y_pan[0] = min(flip0, flip1)
        y_pan[1] = max(flip0, flip1)
    elif x_cart[0] > center[1]:  # Panel 1
        if y_cart[0] < ff1_pix[2]: 
            y_cart[0] = ff1_pix[2]
        if y_cart[1] > ff1_pix[3]: 
            y_cart[1] = ff1_pix[3]
        x_pan = x_cart - ff1_pix[0]
        y_pan = y_cart - ff1_pix[2]

        # Account for horizontal flipping
        midLine = (dexShape[1] - 1) / 2
        flip0 = x_pan[0] + 2 * (midLine - x_pan[0])
        flip1 = x_pan[1] + 2 * (midLine - x_pan[1])
        x_pan[0] = min(flip0, flip1)
        x_pan[1] = max(flip0, flip1)

    # Calculate and return the shape of the ROI
    return y_pan[1] - y_pan[0], x_pan[1] - x_pan[0]

def getEigerPixels(x_cart, y_cart, ff1_pix, eigShape=EIG_SHAPE):
    """
    Computes the pixel coordinates of the ROI for the Eiger detector.

    Parameters:
    -----------
    x_cart : list or ndarray
        X-coordinates of the ROI in the Cartesian system.
    y_cart : list or ndarray
        Y-coordinates of the ROI in the Cartesian system.
    ff1_pix : list or ndarray
        Flat-field parameters for the panel.
    eigShape : tuple, optional
        Shape of the Eiger detector panel (default is (4362, 4148)).

    Returns:
    --------
    tuple
        X and Y coordinates of the ROI in the panel's coordinate system.
    """
    # Ensure coordinates are within bounds
    if x_cart[0] < ff1_pix[0]: 
        x_cart[0] = ff1_pix[0]
    if x_cart[1] > ff1_pix[1]: 
        x_cart[1] = ff1_pix[1]
    if y_cart[0] < ff1_pix[2]: 
        y_cart[0] = ff1_pix[2]
    if y_cart[1] > ff1_pix[3]: 
        y_cart[1] = ff1_pix[3]

    # Transform to panel-specific coordinates
    x_pan = x_cart - ff1_pix[0]
    y_pan = y_cart - ff1_pix[2]
    return x_pan, y_pan

def bilinearInterpMatrix(roiShape, rad_dom, eta_dom, center, detectDist):
    """
    Constructs the bilinear interpolation matrix for mapping Cartesian pixels to polar coordinates.

    Parameters:
    -----------
    roiShape : tuple
        Shape of the ROI (rows, columns).
    rad_dom : ndarray
        Radial domain values.
    eta_dom : ndarray
        Azimuthal domain values.
    center : tuple
        Image center coordinates.
    detectDist : float
        Detector distance.

    Returns:
    --------
    scipy.sparse.coo_array
        Sparse matrix for bilinear interpolation.
    """
    Ainterp = np.zeros((len(rad_dom) * len(eta_dom), roiShape[0] * roiShape[1]))
    k = 0  # Row index in the interpolation matrix

    # Loop through radial and azimuthal domains
    for r in rad_dom:
        for eta in eta_dom:
            # Convert polar to Cartesian coordinates
            x = r * np.cos(eta) + center[1]
            y = -r * np.sin(eta) + center[0]

            # Determine surrounding pixel coordinates
            x1 = np.floor(x)
            x2 = np.ceil(x)
            y1 = np.floor(y)
            y2 = np.ceil(y)

            # Populate the interpolation matrix with weights
            Ainterp[k, int(x2 + y2 * roiShape[1])] = (x - x1) * (y - y1)
            Ainterp[k, int(x2 + y1 * roiShape[1])] = (x2 - x) * (y - y1)
            Ainterp[k, int(x1 + y2 * roiShape[1])] = (x - x1) * (y2 - y)
            Ainterp[k, int(x1 + y1 * roiShape[1])] = (x2 - x) * (y2 - y)
            k += 1

    # Convert to sparse matrix format
    Ainterp = sp.sparse.coo_array(Ainterp)
    return Ainterp

def panelPixelsDex(ff_trans, mmPerPixel, imSize=(4888, 7300), dexShape=DEX_SHAPE):
    """
    Computes pixel boundaries for Dexela detector panels based on flat-field transformations.

    Parameters:
    -----------
    ff_trans : list or ndarray
        Flat-field translation parameters for each panel.
    mmPerPixel : float
        Millimeters per pixel.
    imSize : tuple, optional
        Size of the full detector image (default is (4888, 7300)).
    dexShape : tuple, optional
        Shape of the Dexela detector panel (default is (3888, 3072)).

    Returns:
    --------
    tuple
        Pixel boundaries for panels 1 and 2.
    """
    center = (imSize[0] / 2, imSize[1] / 2)

    # Calculate pixel shifts for each panel
    ff1_xshift = int(round(ff_trans[0] / mmPerPixel))
    ff1_yshift = int(round(ff_trans[1] / mmPerPixel))
    ff2_xshift = int(round(ff_trans[2] / mmPerPixel))
    ff2_yshift = int(round(ff_trans[3] / mmPerPixel))

    # Compute panel boundaries
    ff1_pixels = [
        int(center[1] - dexShape[1] / 2 + ff1_xshift),  # Column start
        int(center[1] + dexShape[1] / 2 + ff1_xshift),  # Column end
        int(center[0] - dexShape[0] / 2 - ff1_yshift),  # Row start
        int(center[0] + dexShape[0] / 2 - ff1_yshift),  # Row end
    ]
    ff2_pixels = [
        int(center[1] - dexShape[1] / 2 + ff2_xshift),  # Column start
        int(center[1] + dexShape[1] / 2 + ff2_xshift),  # Column end
        int(center[0] - dexShape[0] / 2 - ff2_yshift),  # Row start
        int(center[0] + dexShape[0] / 2 - ff2_yshift),  # Row end
    ]

    return ff1_pixels, ff2_pixels

def panelPixelsEiger(ff_trans, mmPerPixel, imSize=(5000, 5000), detectShape=EIG_SHAPE):
    """
    Computes pixel ranges for an Eiger detector panel after translation shifts.

    Parameters:
    -----------
    ff_trans : list or tuple
        Translation shifts [x_shift, y_shift] in millimeters.
    mmPerPixel : float
        Conversion factor from millimeters to pixels.
    imSize : tuple, optional
        Size of the full detector image (height, width).
    detectShape : tuple, optional
        Shape of the detector panel (height, width).

    Returns:
    --------
    ff1_pixels : list
        Pixel boundaries [col_start, col_end, row_start, row_end] for the panel.
    """
    # Determine the center of the image in pixel coordinates
    center = (imSize[0] / 2, imSize[1] / 2)

    # Translate shifts from millimeters to pixels
    ff1_tx = ff_trans[0]
    ff1_ty = ff_trans[1]
    ff1_xshift = int(round(ff1_tx / mmPerPixel))
    ff1_yshift = int(round(ff1_ty / mmPerPixel))

    # Calculate row and column boundaries based on the panel shape and translation
    ff1r1 = int(center[0] - detectShape[0] / 2 - ff1_yshift)
    ff1r2 = int(center[0] + detectShape[0] / 2 - ff1_yshift)
    ff1c1 = int(center[1] - detectShape[1] / 2 + ff1_xshift)
    ff1c2 = int(center[1] + detectShape[1] / 2 + ff1_xshift)

    # Return the calculated pixel boundaries
    ff1_pixels = [ff1c1, ff1c2, ff1r1, ff1r2]
    return ff1_pixels

def polarDomain(detectDist, mmPerPixel, tth, eta, roi_size):
    """
    Constructs radial and azimuthal domains in polar coordinates.

    Parameters:
    -----------
    detectDist : float
        Detector distance from the sample in millimeters.
    mmPerPixel : float
        Conversion factor from millimeters to pixels.
    tth : float
        Two-theta angle in radians.
    eta : float
        Eta angle in radians.
    roi_size : tuple
        Size of the ROI (radial size, angular size).

    Returns:
    --------
    rad_domain : ndarray
        Array of radial positions in pixels.
    eta_domain : ndarray
        Array of azimuthal angles in radians.
    """
    # Calculate the central radial distance in pixels based on two-theta
    r = detectDist * np.tan(tth) / mmPerPixel

    # Define the radial domain (range) based on ROI size
    r1 = r - (roi_size[0] - 1) / 2
    r2 = r + (roi_size[0] - 1) / 2
    rad_domain = np.arange(r1, r2 + 1, 1)

    # Define the azimuthal domain based on ROI size and radial limits
    deta = 1 / r2  # Small angle approximation
    eta1 = eta - roi_size[1] / 2 * deta
    eta2 = eta + roi_size[1] / 2 * deta
    eta_domain = np.linspace(eta1, eta2, roi_size[1])

    return rad_domain, eta_domain

def fetchCartesian(rad_dom, eta_dom, center):
    """
    Converts polar domain boundaries into Cartesian bounds for a given center.

    Parameters:
    -----------
    rad_dom : ndarray
        Array of radial positions.
    eta_dom : ndarray
        Array of azimuthal angles in radians.
    center : tuple
        Center of the image (row, column) in pixel coordinates.

    Returns:
    --------
    x_cart : ndarray
        Minimum and maximum x-coordinates in pixels.
    y_cart : ndarray
        Minimum and maximum y-coordinates in pixels.
    """
    # Compute x-coordinates for both radial boundaries at all eta values
    rad1 = rad_dom[0] * np.cos(eta_dom)
    rad2 = rad_dom[-1] * np.cos(eta_dom)
    x_min1 = np.min(rad1 + center[1])
    x_max1 = np.max(rad2 + center[1])
    x_max2 = np.max(rad1 + center[1])
    x_min2 = np.min(rad2 + center[1])

    # Find overall x-coordinate bounds
    x_min = int(np.floor(min(x_min1, x_min2)))
    x_max = int(np.ceil(max(x_max1, x_max2)))

    # Compute y-coordinates for both radial boundaries at all eta values
    rad1 = rad_dom[0] * np.sin(eta_dom)
    rad2 = rad_dom[-1] * np.sin(eta_dom)
    y_min1 = np.min(-rad1 + center[0])
    y_max1 = np.max(-rad2 + center[0])
    y_max2 = np.max(-rad1 + center[0])
    y_min2 = np.min(-rad2 + center[0])

    # Find overall y-coordinate bounds
    y_min = int(np.floor(min(y_min1, y_min2)))
    y_max = int(np.ceil(max(y_max1, y_max2)))

    return np.array([x_min, x_max + 1]), np.array([y_min, y_max + 1])

def applyTilt(x, y, tilt, detectDist):
    """
    Applies a tilt to Cartesian coordinates.

    Parameters:
    -----------
    x : float
        x-coordinate in pixels.
    y : float
        y-coordinate in pixels.
    tilt : ndarray
        Tilt vector (rotation angles in radians).
    detectDist : float
        Detector distance in millimeters.

    Returns:
    --------
    tuple
        Adjusted (x, y) coordinates after tilt application.
    """
    z = -detectDist  # Set z-coordinate to negative detector distance

    # Create the detector rotation matrix from tilt angles
    R = transforms.xf.makeDetectorRotMat(tilt)

    # Apply rotation matrix to the (x, y, z) vector
    xyz2 = np.array([[x], [y], [z]])
    xyz = np.matmul(R, xyz2)

    return xyz.ravel()[0], xyz.ravel()[1]

def applyTransTilt(x, y, tilt, trans):
    """
    Applies both translation and tilt to Cartesian coordinates.

    Parameters:
    -----------
    x : float
        x-coordinate in pixels.
    y : float
        y-coordinate in pixels.
    tilt : ndarray
        Tilt vector (rotation angles in radians).
    trans : ndarray
        Translation vector (x, y, z) in millimeters.

    Returns:
    --------
    tuple
        Adjusted (x, y) coordinates after applying both translation and tilt.
    """
    z = trans[2]  # Extract z-translation
    x = x - trans[0]  # Subtract x-translation
    y = y - trans[1]  # Subtract y-translation

    # Create the detector rotation matrix from tilt angles
    R = transforms.xf.makeDetectorRotMat(tilt)

    # Apply rotation matrix to the (x, y, z) vector
    xyz2 = np.array([[x], [y], [z]])
    xyz = np.matmul(R, xyz2)

    return xyz.ravel()[0], xyz.ravel()[1]