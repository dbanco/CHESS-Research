# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:21:37 2024

@author: dpqb1
"""
import spotfetch as sf
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import numpy as np
import pickle
import threading
import time
import os

class DataPlotter:
    def __init__(self, root,read_path,spotInds,plotType,titleStr,spotData):
        self.root = root
        self.root.title(titleStr)
        self.read_path = read_path
        self.spotInds = spotInds
        self.plotType = plotType
        self.T = 999999
        # Create figure
        self.fig = Figure(figsize=(10, 5))
        self.gs = gridspec.GridSpec(3, 2, height_ratios=[1,1, 1])
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=5)

        # Create subplots within the same figure
        self.axes = [self.fig.add_subplot(self.gs[i]) for i in range(6)]

        self.trackData = []
        self.dataArray = np.array((6,len(spotInds)))
        self.i = 1
        self.ylabels = [r"$FWHM_\omega$ (deg)",r"$\mu_\omega$ (deg)",\
                        r"$FWHM_\eta$ (deg)",r"$\mu_\eta$ (deg)",\
                        r"$FWHM_{2 \theta}$ (deg)",r"$\mu_{2 \theta}$ (deg)"]

        # self.fig.text(0.5, 0.97, titleStr, ha='center', fontsize=20)
        self.update_plots()
    
    def update_plot(self, event=None):
        # Length of last track
        T = len(self.trackData[-1])
        self.dataArray = np.zeros((6,len(self.spotInds),T))
        self.dataArray[:] = np.nan
        scan = np.zeros((T))
        notDone = True
        # Update all features
        for k in range(len(self.spotInds)):
            if k > len(self.trackData)-1: continue
            for t in range(T):
                if t > len(self.trackData[k])-1: continue
                if len(self.trackData[k][t]) > 0:           
                    avgFWHMeta,avgFWHMtth,avgEta,avgTth = sf.compAvgParams(self.trackData[k][t])
                    FWHMome = sf.estFWHMomega(self.trackData[k][t])
                    Ome = sf.estMEANomega(self.trackData[k][t])
                    self.dataArray[0,k,t] = FWHMome
                    self.dataArray[2,k,t] = avgFWHMeta
                    self.dataArray[4,k,t] = avgFWHMtth
                    self.dataArray[1,k,t] = Ome
                    self.dataArray[3,k,t] = avgEta
                    self.dataArray[5,k,t] = avgTth 
                    # print(f'Spot {k}, Time{t}')

                    if notDone: 
                        scan[t] = self.trackData[k][t][0]['scan']
                        if t == T-1: notDone = False
        
        # Plot each of the features
        for i, ax in enumerate(self.axes):
            ax.clear() 
            if self.plotType == 'Delta':
                for k in range(len(self.spotInds)):
                    ax.plot(scan,self.dataArray[i,k,:]-self.dataArray[i,k,0],'-o')
                ax.set_ylabel(r'$\Delta$ ' + self.ylabels[i])  
            elif self.plotType == 'Mean':
                for k in range(len(self.spotInds)):
                    # print(k)
                    self.dataArray[i,k,:] = sf.mapDiff(self.dataArray[i,k,:]-self.dataArray[i,k,0])
                avg = np.nanmean(self.dataArray[i,:,:],0)
                std = np.nanstd(self.dataArray[i,:,:],0)
                ax.errorbar(scan,avg,std,fmt='-x')
                ax.set_ylabel('Mean ' + self.ylabels[i])  
            else:
                # Orginal
                for k in range(len(self.spotInds)):
                    ax.plot(scan,self.dataArray[i,k,:],'-o')
                ax.set_ylabel(self.ylabels[i])  
            
            ax.set_xlabel("Scan #")  
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        
        self.fig.tight_layout(pad=2.0)
        self.canvas.draw()
            
    def plot(self):
        self.update_plot()    
    
    def update_plots(self):
        def read_data():
            while True:
                try:
                    self.trackData = []
                    self.T = 999999
                    for k in self.spotInds:
                        with open(os.path.join(self.read_path,f'trackData_{k}.pkl'), 'rb') as f:
                            track = pickle.load(f)
                            self.trackData.append(track)
                            self.T = min(self.T,len(self.trackData[-1]))
                    # self.update_plot()
                    self.root.after(0, self.update_plot)  # Schedule update_plot to run in the main thread
                except FileNotFoundError:
                    self.trackData = []
                except pickle.UnpicklingError:
                    self.trackData = []
                except EOFError:
                    self.trackData = []
                time.sleep(1)
        
        threading.Thread(target=read_data, daemon=True).start()
        
def start_gui(read_path,spotInds,plotType,titleStr,spotData):
    root = tk.Tk()
    app = DataPlotter(root,read_path,spotInds,plotType,titleStr,spotData)
    root.mainloop()

if __name__ == "__main__":
    
    topPath = "/mnt/scratch/dbanco/c103_processing/Sample-1/"
    read_path = os.path.join(topPath,'outputs')
    spotData = np.load(os.path.join(topPath,'spots','spots.npz'))

    grains = [276,288,342]
    spotInds = []
    for grain in grains:
        spotInds.append(sf.findSpots(spotData,grains=[grain]))
    titleStr = f'Grain {grains}'
    start_gui(read_path,spotInds,'Mean',titleStr,spotData)