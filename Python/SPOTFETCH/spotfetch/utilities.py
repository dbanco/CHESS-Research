# -*- coding: utf-8 -*-
"""
utilities.py

Description:
Module containing basic functions used across other modules

Functions:
    omegaToFrame
    frameToOmega
    mapOmega
    mapDiff
    wrapFrame
    pathToFile
    timeToFile
    getDataFileSequence
    findSpots
    collectSpotsData
    getSpotID
    matchSpotID
    estMEANomega
    estFWHMomega
    xyToEtaTth
    xyToEtaTthRecenter
    etaTthToPix
    pixToEtaTth

Created on: Fri Apr 12 13:22:47 2024
Author: Daniel Banco
Email: dpqb10@gmail.com
Version: 0.1.0

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
import pandas as pd
import os
import glob
from .detectors import loadYamlData, polarDomain

FRAME1 = 2
NUMFRAMES = 1440
OMEG_RANGE = 360

def omegaToFrame(omega,startFrame=FRAME1,numFrames=NUMFRAMES,omegaRange=OMEG_RANGE,startOmeg = 0):
    '''
    Converts an omega value in radians to its corresponding frame number.

    Parameters:
    -----------
    omega: float
        Omega value in radians, typically in the range (-π, π).
    startFrame: int, optional
        First frame of the scan (default is FRAME1).
    numFrames: int, optional
        Total number of omega frames (default is NUMFRAMES).
    omegaRange: float, optional
        Length of the omega scan in degrees (default is OMEG_RANGE).
    startOmeg: float, optional
        Omega value in radians at `startFrame` (default is 0).

    Returns:
    --------
    frame: int
        Frame number corresponding to the given omega value, ranging 
        from `startFrame` to `startFrame + numFrames`.
    '''
    step = omegaRange/numFrames
    frame = (np.floor((omega*180/np.pi-startOmeg)/step) + startFrame).astype(int)
    return frame
      
def frameToOmega(frame,startFrame=FRAME1,numFrames=NUMFRAMES,omegaRange=OMEG_RANGE,startOmeg = 0):
    '''
    Converts a frame number to its corresponding omega value in radians.
      
    Parameters:
    -----------
    frame: int
        Frame number to be converted.
    startFrame: int, optional
        First frame of the scan (default is FRAME1).
    numFrames: int, optional
        Total number of omega frames (default is NUMFRAMES).
    omegaRange: float, optional
        Length of the omega scan in degrees (default is OMEG_RANGE).
    startOmeg: float, optional
        Omega value in radians at `startFrame` (default is 0).
      
    Returns:
    --------
    omega: float
        Omega value in radians corresponding to the given frame.
    '''
    step = omegaRange/numFrames
    omega = (frame - startFrame)*step + startOmeg
    return omega

def mapOmega(omega):
    '''
    Maps an omega value to the range (-180, 180] degrees.

    Parameters:
    -----------
    omega: float
        Omega value in degrees to be mapped.

    Returns:
    --------
    omega: float
        Omega value mapped to the range (-180, 180].
    '''
    if omega > 180:
        omega = omega - 360
    return omega

def mapDiff(diff):
    '''
    Maps an angular difference to the range (-180, 180] degrees.

    Parameters:
    -----------
    diff: float
        Angular difference in degrees.

    Returns:
    --------
    diff: float
        Angular difference mapped to the range (-180, 180].
    '''
    return (diff + 180) % 360 - 180

def wrapFrame(frm, frm0=FRAME1, numFrms=NUMFRAMES):
    '''
    Wraps a frame number to ensure it stays within the valid frame range.

    Parameters:
    -----------
    frm: int
        Frame number to be wrapped.
    frm0: int, optional
        The starting frame of the scan (default is FRAME1).
    numFrms: int, optional
        Total number of frames in the scan (default is NUMFRAMES).

    Returns:
    --------
    wrappedFrame: int
        Frame number wrapped within the range [frm0, frm0 + numFrms).
    '''
    return np.mod(frm - frm0, numFrms) + frm0

def pathToFile(path):
    '''
    Retrieves and sorts the full file paths of all files in a given directory.

    Parameters:
    -----------
    path: str
        Path to the directory containing the files.

    Returns:
    --------
    fnames: list of str
        A sorted list of full file paths for all files in the specified directory.
    '''
    fnames = os.listdir(path)
    fnames.sort()
    for i in range(len(fnames)):
        fnames[i] = os.path.join(path, fnames[i])
    return fnames

def timeToFile(t, fDir):
    '''
    Constructs file paths for a specific time step by navigating the directory structure.

    Parameters:
    -----------
    t: int
        Time step or identifier for the desired file sequence.
    fDir: str
        Base directory path where time-step-specific directories are located.

    Returns:
    --------
    fnames: list of str
        A sorted list of full file paths for files corresponding to the given time step.
    '''
    dNum = str(t)
    topDir = os.path.join(fDir + dNum, 'ff')
    fnames = pathToFile(topDir)
    return fnames

def getDataFileSequence(dataFile, scanRange):
    '''
    Generates a sequence of data file paths based on a file template and scan range.

    Parameters:
    -----------
    dataFile: str
        Template for the file path, where placeholders can be replaced with values from `scanRange`.
    scanRange: iterable
        A range or list of values to populate the file template with.

    Returns:
    --------
    dataFileSequence: list of str
        A list of full file paths generated from the template and scan range.

    Notes:
    ------
    - The function uses `glob.glob` to find files matching the generated patterns.
    - Assumes each pattern matches exactly one file.
    '''
    dataFileSequence = []
    for scan in scanRange:
        template = dataFile.format(num2=scan)
        pattern = os.path.join(template)
        fname = glob.glob(pattern)[0]
        if len(fname) > 0:
            dataFileSequence.append(fname)
    return dataFileSequence

def findSpots(spotData, grains=None, tths=None, dtth=None, eta=None, deta=None, frm=None):
    '''
    Filters spots from `spotData` based on given conditions such as grain numbers, 
    2theta (tth) values, eta values, and frame indices.

    Parameters:
    -----------
    spotData: dict
        A dictionary containing spot information with the following keys:
        - 'grain_nums': Array of grain numbers.
        - 'etas': Array of eta values.
        - 'tths': Array of two-theta values.
        - 'ome_idxs': Array of frame indices.
    grains: list of int, optional
        List of grain numbers to filter on. Defaults to None (no filtering by grain numbers).
    tths: list of float, optional
        List of two-theta values to filter on. Defaults to None (no filtering by tths).
    dtth: float, optional
        Tolerance for filtering by tths. Required if `tths` is provided.
    eta: float, optional
        Central eta value for filtering. Defaults to None (no filtering by eta).
    deta: float, optional
        Tolerance for filtering by eta. Required if `eta` is provided.
    frm: int, optional
        Frame index to filter on. Defaults to None (no filtering by frame indices).

    Returns:
    --------
    spotInds: ndarray
        Array of indices corresponding to the spots that satisfy the conditions.
    '''
    grain_nums = spotData.get('grain_nums', None)
    etas = spotData.get('etas', None)
    tths_data = spotData.get('tths', None)
    frms = spotData.get('ome_idxs', None)
    
    # Conditions
    cond1 = np.isin(grain_nums, grains) if grains is not None else True
    cond2 = ((etas > eta - deta) & (etas < eta + deta)) if eta is not None else True
    cond3 = (
        np.any([(tths_data > tth - dtth) & (tths_data < tth + dtth) for tth in tths], axis=0) 
        if tths is not None and dtth is not None else True
    )
    cond4 = (frms == frm) if frm is not None else True

    # Combine all conditions
    spotInds = np.where(cond1 & cond2 & cond3 & cond4)[0]

    return spotInds

def collectSpotsData(outPath, spotsPath, sortFlag=False):
    '''
    Collects and processes spot data from multiple `.out` files in a given directory structure, 
    aggregates the data, filters invalid entries, and saves the result as a `.npz` file.

    Parameters:
    -----------
    outPath: str
        Path to the directory where the output `.npz` file will be saved.
    spotsPath: str
        Path to the directory containing the subdirectories with `.out` files.
    sortFlag: bool, optional
        If True, sorts the data by omega values before saving. Default is False.

    Process:
    --------
    - Reads `.out` files within subdirectories of `spotsPath`.
    - Extracts relevant data fields from each file (e.g., coordinates, grain numbers, angles).
    - Appends data across files and filters invalid entries (e.g., missing or NaN values).
    - Optionally sorts the data by omega values (`omes`) if `sortFlag` is set.
    - Saves the processed data to a compressed `.npz` file.

    Output:
    -------
    Saves a `.npz` file in `outPath` containing:
    - `Xs`, `Ys`: Predicted X and Y coordinates.
    - `Xm`, `Ym`: Measured X and Y coordinates.
    - `id_nums`: IDs of the spots.
    - `tths`, `etas`, `omes`: Measured two-theta, eta, and omega values.
    - `tths_pred`, `etas_pred`, `omes_pred`: Predicted two-theta, eta, and omega values.
    - `ome_idxs`: Frame indices corresponding to omega values.
    - `grain_nums`: Grain numbers for each spot.
    - `PID`: Phase IDs.
    - `H`, `K`, `L`: Miller indices.

    Notes:
    ------
    - Files with invalid or missing data are automatically filtered out.
    - Assumes `.out` files have specific columns such as 'pred X', 'meas X', 'meas tth', etc.
    '''
    all_entries = os.listdir(spotsPath)
    directories = [entry for entry in all_entries if os.path.isdir(os.path.join(spotsPath, entry))]
    fold_names = sorted(directories)

    created = False

    for fold_name in fold_names:
        file_names = sorted(os.listdir(os.path.join(spotsPath, fold_name)))
        print(fold_name)

        for file_name in file_names:
            if file_name.endswith(".out"):
                file_path = os.path.join(spotsPath, fold_name, file_name)

                # Load .out file
                df = pd.read_csv(file_path, sep='\\s+', engine='python')

                # Extract data
                grain_number = int(file_name[-7:-4])
                new_data = {
                    'Xs': df['pred X'].to_numpy(),
                    'Ys': df['pred Y'].to_numpy(),
                    'Xm': df['meas X'].to_numpy(),
                    'Ym': df['meas Y'].to_numpy(),
                    'id_nums': df['# ID'].to_numpy(),
                    'PID': df['PID'].to_numpy(),
                    'H': df['H'].to_numpy(),
                    'K': df['K'].to_numpy(),
                    'L': df['L'].to_numpy(),
                    'grain_nums': grain_number * np.ones(df.shape[0]),
                    'tths': df['meas tth'].to_numpy(),
                    'etas': df['meas eta'].to_numpy(),
                    'omes': df['meas ome'].to_numpy(),
                    'tths_pred': df['pred tth'].to_numpy(),
                    'etas_pred': df['pred eta'].to_numpy(),
                    'omes_pred': df['pred ome'].to_numpy(),
                }

                # Initialize or append
                if not created:
                    combined_data = new_data
                    created = True
                else:
                    for key in combined_data:
                        combined_data[key] = np.append(combined_data[key], new_data[key])

    # Filter invalid data
    invalid_mask = (
        (combined_data['id_nums'] == -999) | 
        np.isnan(combined_data['tths']) | 
        np.isnan(combined_data['etas'])
    )
    for key in combined_data:
        combined_data[key] = np.delete(combined_data[key], invalid_mask)

    # Convert omegas to frame indices
    combined_data['ome_idxs'] = omegaToFrame(combined_data['omes'])

    # Sort data if required
    if sortFlag:
        sort_indices = np.argsort(combined_data['omes'])
        for key in combined_data:
            combined_data[key] = combined_data[key][sort_indices]

    # Save to .npz file
    save_file = os.path.join(outPath, 'spots.npz')
    np.savez(save_file, **combined_data)

def getSpotID(spotData, k):
    """
    Retrieves the spot ID for a specific index `k` from the spot data.

    Parameters:
    -----------
    spotData: dict
        Dictionary containing spot data with keys like 'grain_nums', 'PID', etc.
    k: int
        Index of the spot.

    Returns:
    --------
    list
        List containing [grain number, PID, H, K, L, ID].
    """
    keys = ['grain_nums', 'PID', 'H', 'K', 'L', 'id_nums']
    spot_id = [spotData[key][k] for key in keys]
    return spot_id

def matchSpotID(spotData, spot_id, k):
    """
    Finds the index of a spot in `spotData` that matches the given `spot_id`.

    Parameters:
    -----------
    spotData: dict
        Dictionary containing spot data with keys like 'grain_nums', 'PID', etc.
    spot_id: list
        List containing [grain number, PID, H, K, L, ID] of the spot to match.
    k: int
        Original index of the spot for reference in case of no match.

    Returns:
    --------
    int or float
        Index of the matching spot, or NaN if no match is found.
    """
    # Extract spot data fields
    fields = ['grain_nums', 'PID', 'H', 'K', 'L', 'id_nums']
    grNum, PID, H, K, L, idNum = [spotData[field] for field in fields]

    # Match based on grain number, PID, and Miller indices
    bin_array = (grNum == spot_id[0]) & (PID == spot_id[1]) & \
                (H == spot_id[2]) & (K == spot_id[3]) & (L == spot_id[4])
    match_ind = np.where(bin_array)[0]

    # Refine match if multiple indices found
    if len(match_ind) > 1:
        diffs = abs(idNum[match_ind] - spot_id[5])
        valid_matches = diffs < 50
        match_ind = match_ind[valid_matches]
        diffs = diffs[valid_matches]

        if len(match_ind) > 1:
            match_ind = match_ind[np.argmin(diffs)]

    # Handle no matches
    if not match_ind.size:
        print(f'Match not found for spot {k}')
        return np.nan

    return match_ind.item()  # Return scalar value

def estMEANomega(track):
    """
    Estimates the mean omega value for a given track.

    Parameters:
    -----------
    track: list of dicts
        Each dictionary contains keys 'roi' (region of interest array) 
        and 'frm' (frame index).

    Returns:
    --------
    float
        Estimated mean omega value.
    """
    step = OMEG_RANGE / NUMFRAMES

    if len(track) == 1:
        return frameToOmega(track[0]['frm'])

    roiOmega = np.array([np.nansum(t['roi']) for t in track])
    omegaRange = np.array([frameToOmega(t['frm']) for t in track])
    weighted_mean_index = np.sum(np.arange(len(track)) * roiOmega) / np.sum(roiOmega)

    ind1 = int(np.floor(weighted_mean_index))
    ind2 = int(np.ceil(weighted_mean_index))

    meanOmega = (
        omegaRange[ind1] + (weighted_mean_index - ind1) * step
        if omegaRange[ind1] > omegaRange[ind2]
        else weighted_mean_index * step
    )

    return mapOmega(meanOmega)  # Map to -180° to 0° if necessary

def estFWHMomega(track):
    """
    Estimates the Full-Width Half-Maximum (FWHM) of omega for a given track.

    Parameters:
    -----------
    track: list of dicts
        Each dictionary contains keys 'roi' (region of interest array) 
        and 'frm' (frame index).

    Returns:
    --------
    float
        Estimated FWHM of omega.
    """
    step = OMEG_RANGE / NUMFRAMES

    if len(track) == 1:
        return 0

    roiOmega = np.array([np.nansum(t['roi']) for t in track])
    weighted_mean_index = np.sum(np.arange(len(track)) * roiOmega) / np.sum(roiOmega)

    # Compute variance
    varOmega = np.sum(
        roiOmega * ((np.arange(len(track)) - weighted_mean_index) ** 2)
    ) / np.sum(roiOmega) * step**2

    # Calculate FWHM
    fwhmOmega = 2 * np.sqrt(2 * np.log(2) * varOmega)

    return fwhmOmega

def xyToEtaTth(x, y, params):
    """
    Converts Cartesian coordinates (x, y) to eta and 2-theta.

    Parameters:
    -----------
    x, y: float or np.ndarray
        Cartesian coordinates.
    params: dict
        Dictionary containing detector parameters.

    Returns:
    --------
    tuple
        eta (float or np.ndarray): Azimuthal angle in radians.
        tth (float or np.ndarray): 2-theta angle in radians.
    """
    detectDist, _, _, _ = loadYamlData(params)
    eta = np.arctan2(y, x)
    rad = np.sqrt(x**2 + y**2)
    tth = np.arctan(rad / detectDist)
    return eta, tth

def xyToEtaTthRecenter(x, y, params):
    """
    Converts Cartesian coordinates (x, y) to eta and 2-theta after recentering.

    Parameters:
    -----------
    x, y: float or np.ndarray
        Cartesian coordinates.
    params: dict
        Dictionary containing detector parameters, including translation.

    Returns:
    --------
    tuple
        eta (float or np.ndarray): Azimuthal angle in radians.
        tth (float or np.ndarray): 2-theta angle in radians.
    """
    detectDist, _, trans, _ = loadYamlData(params)
    x += trans[0]
    y += trans[1]
    eta = np.arctan2(y, x)
    rad = np.sqrt(x**2 + y**2)
    tth = np.arctan(rad / detectDist)
    return eta, tth
    
def etaTthToPix(eta, tth, etaRoi, tthRoi, params):
    """
    Converts eta and 2-theta to pixel coordinates in the detector.

    Parameters:
    -----------
    eta, tth: float or np.ndarray
        Eta (azimuthal) and 2-theta angles in radians.
    etaRoi, tthRoi: tuple
        Region of interest for eta and 2-theta, respectively.
    params: dict
        Dictionary containing detector parameters, including ROI size.

    Returns:
    --------
    tuple
        row_pos, col_pos: Pixel coordinates.
    """
    roiSize = params['roiSize']
    detectDist, mmPerPixel, _, _ = loadYamlData(params, tthRoi, etaRoi)
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tthRoi, etaRoi, roiSize)

    # Calculate radial and azimuthal positions
    rad = np.tan(tth) * detectDist / mmPerPixel
    row_pos = (rad - rad_dom[0]) / (rad_dom[-1] - rad_dom[0]) * (roiSize[0] - 1)
    col_pos = (eta - eta_dom[0]) / (eta_dom[-1] - eta_dom[0]) * (roiSize[1] - 1)

    return row_pos, col_pos
    
def pixToEtaTth(p1, p2, tthRoi, etaRoi, params):
    """
    Converts pixel coordinates to eta and 2-theta.

    Parameters:
    -----------
    p1, p2: float
        Pixel coordinates.
    tthRoi, etaRoi: tuple
        Region of interest for 2-theta and eta, respectively.
    params: dict
        Dictionary containing detector parameters, including ROI size.

    Returns:
    --------
    tuple
        etaNew, tthNew: New eta and 2-theta values.
        deta, dtth: Step sizes in eta and 2-theta.
    """
    roiSize = params['roiSize']
    detectDist, mmPerPixel, _, _ = loadYamlData(params, tthRoi, etaRoi)
    rad_dom, eta_dom = polarDomain(detectDist, mmPerPixel, tthRoi, etaRoi, roiSize)

    # Calculate step sizes
    deta = abs(eta_dom[1] - eta_dom[0])
    hypot = detectDist * np.cos(tthRoi)
    dtth = np.arctan(mmPerPixel / hypot)

    # Determine new eta and radial positions
    i1, j1 = int(np.floor(p1)), int(np.floor(p2))
    etaNew = eta_dom[i1] + deta * (p1 % 1)
    radNew = rad_dom[j1] + (p2 % 1)
    tthNew = np.arctan(radNew * mmPerPixel / detectDist)

    return etaNew, tthNew, deta, dtth
     
# def visualTrackData(trackData):
#     T = len(trackData)
#     K = len(trackData[0])
#     FWHMeta = np.zeros((T,K))
#     meaneta = np.zeros((T,K))
    
#     for t in range(T):
#         if trackData[t] != []:
#             for k in range(K):
#                 if trackData[t][k] != []:
#                     FWHMeta[t] = trackData[t][k][0]['p'][3]
#                     meaneta[t] = trackData[t][k][0]['eta']
#     return FWHMeta,meaneta
            

