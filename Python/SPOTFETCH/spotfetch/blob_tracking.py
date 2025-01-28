# -*- coding: utf-8 -*-
"""
spotfetch.blob_tracking

Description:
Tools for tracking spots using 3D blobs in X-ray diffraction data

Functions:
trackSpotBlob
initSpotBlob
processSpotBlob
compareBlobs
assembleBlobTrack
DoG
detectBlobDoG
blobFeaturesDoG

Created on Wed Jan 22 18:49:48 2025
Author: Daniel Banco
Email: dpqb10@gmail.com
Version: 0.1.0
"""
import numpy as np
import pickle
import os
import time
import spotfetch as sf
from scipy.ndimage import gaussian_filter
from scipy.ndimage import center_of_mass
from skimage.feature import hessian_matrix
from scipy.ndimage import label

def trackSpotBlob(spotInd, spotData, dataFileSequence, trackPath, params):
    """
    Tracks a spot blob across a sequence of data files.

    Parameters:
    -----------
    spotInd: int
        Index of the spot in the spot data.
    spotData: dict
        Dictionary containing information about detected spots.
        Must include keys: 'Xm', 'Ym', 'ome_idxs'.
    dataFileSequence: list
        List of file names for the data sequence to process.
    trackPath: str
        Directory path to save tracking data files.
    params: dict
        Dictionary of parameters required for spot tracking.

    Returns:
    --------
    int
        Returns 0 if tracking was skipped because the file already exists.
        Otherwise, returns None.
    """
    # Define the track file path
    trackFile = os.path.join(trackPath, f'trackData_{spotInd}.pkl')

    # Check if the tracking file already exists
    if os.path.exists(trackFile):
        print(f'File exists. Spot {spotInd} not tracked.')
        return 0

    # Initial spot location
    x = spotData['Xm'][spotInd]
    y = spotData['Ym'][spotInd]
    frm = int(spotData['ome_idxs'][spotInd])

    # Convert to eta and 2-theta coordinates
    eta, tth = sf.xyToEtaTthRecenter(x, y, params)

    # Tracking across the data sequence
    for t_ind, fileName in enumerate(dataFileSequence):
        if t_ind == 0:
            # Initialize the spot blob for the first frame
            initSpotBlob(spotInd, eta, tth, frm, fileName, trackFile, params)
        else:
            # Process subsequent frames
            processSpotBlob(spotInd, t_ind, fileName, trackFile, params)

def initSpotBlob(spotInd, eta, tth, frm, fileName, trackFile, params):
    """
    Initializes spot blob tracking by detecting and analyzing a region of interest (ROI).
    Saves tracking outputs to the designated trackFile

    Parameters:
    -----------
    spotInd: int
        Index of the spot to initialize tracking.
    eta: float
        Eta coordinate of the spot.
    tth: float
        2-theta (tth) coordinate of the spot.
    frm: int
        Frame index of the spot.
    fileName: str
        File name containing the data for the spot.
    trackFile: str
        Path to save the tracking data.
    params: dict
        Dictionary containing tracking parameters.

    Returns:
    --------
    None
    """
    if params.get('benchmarkFlag', False):
        start_time = time.time()

    # Initialize output data structure
    outData = {
        'trackData': [],
        'trackTimes': []
    }

    # 1. Load 3D region of interest (ROI)
    roi3D = sf.loadPolarROI3D(fileName, tth, eta, frm, params)

    # 2. Perform blob detection
    blobs, num_blobs, hess_mat = sf.detectBlobDoG(roi3D)

    # 3. Extract blob features
    features = sf.blobFeaturesDoG(roi3D, blobs, num_blobs, hess_mat)

    # 4. Identify the initial spot blob closest to the center of ROI
    center = (np.array(params['roiSize']) - 1) / 2
    distances = np.sum((np.array(features['com']) - center) ** 2, axis=1)
    closest_blob_idx = np.argmin(distances)

    # Visualize the ROI and detected blobs for debugging
    # sf.showTensor(roi3D)
    # sf.showTensor(blobs)

    # Check if the closest blob is within the initialization distance threshold
    if distances[closest_blob_idx] < params.get('initDist', np.inf):
        # Isolate the blob and create a new track
        blob = roi3D.copy()
        blob[~(blobs == (closest_blob_idx + 1))] = 0

        newTrack = assembleBlobTrack(tth, eta, frm, 0, blob,
                                     features[closest_blob_idx], params)

        outData['trackData'].append(newTrack)

        # Visualize the isolated blob for debugging
        sf.showTensor(blob)
    else:
        # If no valid blob is found, append False
        outData['trackData'].append(False)

    # Benchmarking: Record processing time
    if params.get('benchmarkFlag', False):
        end_time = time.time()
        outData['trackTimes'].append(end_time - start_time)

    # Save the output data to a file
    with open(trackFile, 'wb') as f:
        pickle.dump(outData, f)

def processSpotBlob(spotInd, t_ind, fnames, trackFile, params):
    """
    Processes and tracks a spot blob for a given frame. Updates the track data and 
    appends new information to the track file.

    Parameters:
    -----------
    spotInd: int
        Index of the spot being processed.
    t_ind: int
        Index of the current time step in the sequence.
    fnames: list of str
        File names corresponding to the data for each frame.
    trackFile: str
        Path to the file storing track data.
    params: dict
        Dictionary of parameters for processing and tracking.
        
    Returns:
    --------
    int
        Returns 0 if tracking is skipped; otherwise, no return value.
    """
    # Start benchmark timer if enabled
    if params.get('benchmarkFlag', False):
        start_time = time.time()

    # Load existing track data
    with open(trackFile, 'rb') as f:
        outData = pickle.load(f)

    # If track data is empty, exit early
    if len(outData['trackData']) == 1 and outData['trackData'][0] == []:
        return 0

    # Only continue tracking if the spot hasn't been lost
    if outData['trackData'][-1] == False:
        outData['trackData'].append(False)
    else:
        # 1. Load the previous track's data
        prevTrack = outData['trackData'][-1]
        prevTth = prevTrack['tth']
        prevEta = prevTrack['eta']
        prevFrm = prevTrack['frm']

        # 2. Load ROI for the current frame
        roi3D = sf.loadPolarROI3D(fnames, prevTth, prevEta, prevFrm, params)

        # 3. Perform blob detection
        blobs, num_blobs, hess_mat = sf.detectBlobDoG(roi3D)

        # 4. Extract blob features
        features = sf.blobFeaturesDoG(roi3D, blobs, num_blobs, hess_mat)

        # 5. Compare detected blobs to previous track
        selInd = compareBlobs(
            prevTrack, prevTth, prevEta, prevFrm, t_ind, features, params
        )

        # 6. If a matching blob is found, update the track
        if selInd == -1:
            outData['trackData'].append(False)  # Spot lost
        else:
            # Update blob for current track
            blob = roi3D.copy()
            blob[~(blobs == selInd)] = 0
            newTrack = assembleBlobTrack(
                prevTth, prevEta, prevFrm, t_ind, blob,
                features['RT'][selInd],features['ST'][selInd],
                features['AT'][selInd],features['sum'][selInd],
                features['com'][selInd], params
            )
            outData['trackData'].append(newTrack)

    # Benchmarking end and time recording
    if params.get('benchmarkFlag', False):
        end_time = time.time()
        outData['trackTimes'].append(end_time - start_time)

    # Save the updated track data to the track file
    with open(trackFile, 'wb') as f:
        pickle.dump(outData, f)

def compareBlobs(prevTrack, prevTth, prevEta, prevFrm, t_ind, features, params):
    """
    Compares blobs to the previous track and selects the best match.

    Parameters:
    -----------
    prevTrack : dict
        The previous track information.
    prevTth : float
        Previous theta value.
    prevEta : float
        Previous eta value.
    prevFrm : int
        Previous frame index.
    t_ind : int
        Current time index.
    features : dict
        Dictionary containing blob features (RT, ST, AT, sum, com).
    params : dict
        Dictionary of parameters, including thresholds and ROI size.

    Returns:
    --------
    int
        Index of the selected blob or -1 if no match is found.
    """
    com = features['com']
    scores = np.full(len(com), np.inf)  # Initialize scores with infinity
    
    for i, center in enumerate(com):
        if np.isnan(center).any():
            continue
        
        etaNew, tthNew, deta, dtth = sf.pixToEtaTth(center[1], center[0], prevTth, prevEta, params)
        frmNew = prevFrm + center[2] - (params['roiSize'][2] - 1) / 2
        
        # 1. Check if within minimum distance
        tthCheck = abs(tthNew - prevTth) < params['gamma'][0] * dtth
        etaCheck = abs(etaNew - prevEta) < params['gamma'][1] * deta
        frmCheck = abs(frmNew - prevFrm) < params['gamma'][2]
        
        if tthCheck and etaCheck and frmCheck:
            # 2. Compute match score
            RTdiff = abs(prevTrack['RT'] - features['RT'][i])
            STdiff = abs(prevTrack['ST'] - features['ST'][i]) / prevTrack['ST']
            ATdiff = abs(prevTrack['AT'] - features['AT'][i]) / prevTrack['AT']
            sumIntdiff = abs(prevTrack['sum'] - features['sum'][i]) / prevTrack['sum']
            blobDist = np.sqrt(
                ((tthNew - prevTth) / dtth)**2 +
                ((etaNew - prevEta) / deta)**2 +
                (frmNew - prevFrm)**2
            ) / np.mean(params['gamma'][:2])
            
            scores[i] = RTdiff + STdiff + ATdiff + sumIntdiff + blobDist
    
    selInd = np.argmin(scores)
    return -1 if scores[selInd] == np.inf else selInd
        
def assembleBlobTrack(tthRoi, etaRoi, frmRoi, scan, blob, features, params):
    """
    Assembles a new track from the identified blob.

    Parameters:
    -----------
    tthRoi, etaRoi, frmRoi : float
        ROI coordinates (theta, eta, and frame index).
    scan : int
        Current scan index.
    blob : ndarray
        The blob data.
    features : dict
        Blob features.
    params : dict
        Parameters for feature extraction.

    Returns:
    --------
    dict
        New track dictionary with updated information.
    """
    com = features['com']
    etaNew, tthNew, deta, dtth = sf.pixToEtaTth(com[1], com[0], tthRoi, etaRoi, params)
    frmNew = frmRoi + com[2] - (params['roiSize'][2] - 1) / 2
    
    return {
        'etaRoi': etaRoi,
        'tthRoi': tthRoi,
        'frmRoi': frmRoi,
        'eta': etaNew,
        'tth': tthNew,
        'frm': frmNew,
        'scan': scan,
        'deta': deta,
        'dtth': dtth,
        'blob': blob,
        'features': features
    }
          
def DoG(f, sigma, dsigma, gamma=2):
    """
    Computes the Difference of Gaussians (DoG) approximation.

    Parameters:
    -----------
    f : ndarray
        Input image.
    sigma : float
        Base Gaussian sigma.
    dsigma : float
        Increment for sigma.
    gamma : float, optional
        Scaling factor, default is 2.

    Returns:
    --------
    ndarray
        Normalized DoG result.
    """
    g1 = gaussian_filter(f, sigma=sigma)
    g2 = gaussian_filter(f, sigma=sigma + dsigma)
    return (g2 - g1) / (sigma * dsigma)

def detectBlobDoG(x):
    """
    Detects blobs using the Difference of Gaussians (DoG) method.

    Parameters:
    -----------
    x : ndarray
        Input 3D image.

    Returns:
    --------
    tuple
        - blobs : ndarray
            Labeled blob regions.
        - num_blobs : int
            Number of blobs detected.
        - hess_mat : list of ndarray
            Hessian matrices of the DoG result.
    """
    # 1. Compute normalized DoG
    dog_norm = DoG(x, sigma=2, dsigma=1.5)

    # 2. Pre-segmentation
    hess_mat = hessian_matrix(dog_norm)
    D1 = np.zeros(hess_mat[0].shape)
    D2 = np.zeros(hess_mat[0].shape)
    D3 = np.zeros(hess_mat[0].shape)
    for i1 in range(hess_mat[0].shape[0]):
        for i2 in range(hess_mat[0].shape[1]):
            for i3 in range(hess_mat[0].shape[2]):
                h_mat = np.array([
                    [hess_mat[0][i1, i2, i3], hess_mat[1][i1, i2, i3], hess_mat[2][i1, i2, i3]],
                    [hess_mat[1][i1, i2, i3], hess_mat[3][i1, i2, i3], hess_mat[4][i1, i2, i3]],
                    [hess_mat[2][i1, i2, i3], hess_mat[4][i1, i2, i3], hess_mat[5][i1, i2, i3]]
                ])
                D1[i1,i2,i3] = h_mat[0,0]
                D2[i1,i2,i3] = np.linalg.det(h_mat[:2,:2])
                D3[i1,i2,i3] = np.linalg.det(h_mat)

    posDefIndicator = (D1 > 0) & (D2 > 0) & (D3 > 0)
    blobs, num_blobs = label(posDefIndicator)
    return blobs, num_blobs, hess_mat

def blobFeaturesDoG(x, blobs, num_blobs, hess_mat):
    """
    Extracts blob features from the given 3D ROI data.

    Parameters:
    -----------
    x: ndarray
        3D HEXD data (e.g., intensity values).
    blobs: ndarray
        Blob labels (each blob has a unique label).
    num_blobs: int
        Number of blobs detected.
    hess_mat: list of ndarray
        Hessian matrices for each component of the blob.

    Returns:
    --------
    dict
        A dictionary containing the following features:
        - Blobness feature (RT)
        - Flatness feature (ST)
        - Average intensity feature (AT)
        - Total intensity feature (sum)
        - Center of mass feature (com)
    """
    
    # Initialize lists to hold the features for all blobs
    RT, ST, AT, sum_intensity, com = [], [], [], [], []
    
    for i in range(1, num_blobs + 1):
        # Find indices of the current blob
        b_indices = blobs == i
        hess_i = []
        
        # Sum Hessian components for the current blob
        for component in hess_mat:
            hess_i.append(np.sum(component[b_indices]))
        
        # Construct the Hessian matrix for eigenvalue calculation
        h_mat = np.array([[hess_i[0], hess_i[1], hess_i[2]],
                          [hess_i[1], hess_i[3], hess_i[4]],
                          [hess_i[2], hess_i[4], hess_i[5]]])
        
        # Calculate the eigenvalues
        eigs_i = np.linalg.eig(h_mat)
        
        # Blobness feature
        RT_value = 3 * np.abs(eigs_i[0][0] * eigs_i[0][1] * eigs_i[0][2])**(2/3) / (
            np.abs(eigs_i[0][0] * eigs_i[0][1]) + 
            np.abs(eigs_i[0][0] * eigs_i[0][2]) + 
            np.abs(eigs_i[0][2] * eigs_i[0][1])
        )
        
        # Flatness feature
        ST_value = np.sqrt(eigs_i[0][0]**2 + eigs_i[0][1]**2 + eigs_i[0][2]**2)
        
        # Average intensity feature
        AT_value = np.mean(x[b_indices])
        
        # Total intensity feature
        sumIntensity = np.mean(x[b_indices])
        
        # Center of mass feature
        xx = x.copy()
        xx[~b_indices] = 0
        com_value = center_of_mass(xx)
        
        # Append the features to the respective lists
        RT.append(RT_value)
        ST.append(ST_value)
        AT.append(AT_value)
        sum_intensity.append(sumIntensity)
        com.append(com_value)
    
    # Return a dictionary containing the features
    return {
        'RT': RT,
        'ST': ST,
        'AT': AT,
        'sum': sum_intensity,
        'com': com
    }

# def roiTensor(spotInd,spotData,dome,scanRange,dataPath,spotFiles,params):

#     T = len(scanRange)
#     M_ome = 2*dome+1
  
#     # Load sim truth
#     spotDataList = []
#     for i in range(len(spotFiles)):
#         if os.path.exists(spotFiles[i]):
#             spotDataList.append(np.load(spotFiles[i]))
    
#     # Get initial spot location
#     x = spotData['Xm'][spotInd]
#     y = spotData['Ym'][spotInd]
#     eta0, tth0 = sf.xyToEtaTthRecenter(x,y,params)
#     frm0 = spotData['ome_idxs'][spotInd]
#     print(spotData['omes'][spotInd])
    
#     etaRoi = eta0
#     tthRoi = tth0
#     frmRange = np.arange(frm0-dome,frm0+dome+1)

#     M_tth = params['roiSize'][0]
#     M_eta = params['roiSize'][1]
#     roiTensor = np.zeros((M_tth,M_eta,M_ome,T))