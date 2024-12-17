# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 10:17:23 2024

@author: dpqb1
"""
import spotfetch as sf
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.gridspec as gridspec
import numpy as np
import pickle
import threading
import time
import os

read_path = "C:\\Users\\dpqb1\\Documents\\Data\\ti-2\\outputs\\"
with open(os.path.join(read_path,f'trackData_0.pkl'), 'rb') as f:
    trackData = pickle.load(f)
    
spotInds = [0,1,2,3]
K = len(spotInds)
T = len(trackData)
omegTracks = np.zeros((T,K))
avgFWHMeta = np.zeros((T,K))
avgFWHMtth = np.zeros((T,K))
avgEta     = np.zeros((T,K))
avgTth     = np.zeros((T,K))

for k in spotInds:
    # Open spot k track
    with open(os.path.join(read_path,f'trackData_{k}.pkl'), 'rb') as f:
        trackData = pickle.load(f)

    for t in range(T):
        if not (trackData[t] == []):
            omegTracks[t,k] = len(trackData[t])
            fwhm1, fwhm2, aEta, aTth = sf.compAvgParams(trackData[t])
            avgFWHMeta[t,k] = fwhm1
            avgFWHMtth[t,k] = fwhm2
            avgEta[t,k] = aEta
            avgTth[t,k] = aTth
            # print(f'avgFWHMeta = {avgFWHMeta[t,k]}')
            
plt.figure()
plt.subplot(5,1,1)
plt.plot(avgFWHMeta)
plt.subplot(5,1,2)
plt.plot(avgFWHMtth)
plt.subplot(5,1,3)
plt.plot(avgEta)
plt.subplot(5,1,4)
plt.plot(avgTth)
plt.subplot(5,1,5)
plt.plot(omegTracks)
