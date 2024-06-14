# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:45:52 2024

@author: dpqb1
"""
import pickle
import spotfetch as sf
import matplotlib.pyplot as plt

k = 0
outPath = "C:\\Users\\dpqb1\\Documents\\Data\\ti-2\\outputs\\"
track_file = outPath + f'trackData_{k}.pkl'

i = 5
t = 56
dataPath = "D:\\CHESS_raw_data\\ti-2-tension\\"
dataDir = dataPath + f'{t}\\ff\\'

params = {}
# params['detector'] = 'eiger'
params['detector'] = 'dexela'
params['yamlFile'] = "C:\\Users\\dpqb1\\Documents\\Data\\ti-2\\dex-refined-1.yml"
params['imSize'] = (4888,7300) 
params['roiSize'] = 40
params['gamma'] = [4,8,6,6]

with open(track_file, 'rb') as f:
    track_data = pickle.load(f)
    
track_data[i]

prevTracks = track_data[i-1]

fnames = sf.pathToFile(dataDir)
fname1 = dataDir + fnames[0]
fname2 = dataDir + fnames[1]

eta = prevTracks[0]['eta']
tth = prevTracks[0]['tth']
frm = prevTracks[0]['frm']

newTrack, peakFound = sf.evaluateROI(fname1,fname2,prevTracks,\
                    tth,eta,1444,t,params)
    
plt.figure()
plt.subplot(3,1,1)
plt.imshow(prevTracks[0]['roi'])
plt.subplot(3,1,2)
plt.imshow(track_data[i][0]['roi'])
plt.subplot(3,1,3)
plt.imshow(newTrack['roi'])