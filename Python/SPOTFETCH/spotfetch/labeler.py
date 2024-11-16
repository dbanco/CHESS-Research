#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spotfetch.labeler

Utilities for labeling and annotating X-ray diffraction data.

Created on: Fri Nov 15 23:00:28 2024
Author: Daniel Banco
"""
import pickle
import os
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RectangleSelector
from tkinter import font
import detectors as dt

class truthLabeler:
    def __init__(self,root,spotInd,scan_ind,spotData,scanRange,\
                 trackPath,truthPath,dataPath,params):

        self.root = root
        self.root.title("Data Labeler")
        self.root.geometry("1200x1500")
        self.root.resizable(True, True)
        
        # Initialize plot
        plt.ioff()
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.X, expand=True)
        
        self.eta0 = spotData['etas'][spotInd]
        self.tth0 = spotData['tths'][spotInd]
        self.frm0 = spotData['ome_idxs'][spotInd] 
        self.frm = spotData['ome_idxs'][spotInd]  
        self.etaRoi = self.eta0
        self.tthRoi = self.tth0
    
        self.scan_ind = scan_ind
        self.scan = scanRange[self.scan_ind]
        
        self.dataPath = dataPath
        self.scanRange = scanRange
        self.params = params
        self.trackFound = False,
        self.truthFound = False,
        self.om_ind1 = []
        self.om_ind2 = []
        
        helpMenu = 'Controls: \n'+\
              '--------------------\n'+\
              'w: -1 Omega frame \n'+\
              's: +1 Omega frame \n'+\
              'a: -1 Scan (t) \n'+\
              'd: +1 Scan (t) \n'+\
              'x: Finish \n'+\
              'p: Enter truth position \n'+\
              'q: Save truth \n'+\
              '9: Remove current truth \n'+\
              'h: Help menu'
        self.helpMenu = helpMenu
        
        # Load track data
        self.track_file = os.path.join(trackPath,f'trackData_{spotInd}.pkl')
        if os.path.exists(self.track_file):
            with open(self.track_file, 'rb') as f:
                self.track_data = pickle.load(f)
                print('Track loaded')
        else:
            self.track_data = []
            print('No track data')
            
        # Load truth data 
        self.truth_file = os.path.join(trackPath,f'truthData_{spotInd}.pkl')
        if os.path.exists(self.truth_file):
            with open(self.truth_file, 'rb') as f:
                self.truth_data = pickle.load(f)
                print('Truth loaded')
        else:   
            T = len(scanRange)
            self.truth_data = [[] for _ in range(T)]
            print('Creating new truth')
        
        bold_font = font.Font(size=16, weight="bold")
        # Create frames to control widget arrangement
        left_frame = tk.Frame(root)
        right_frame = tk.Frame(root)
        third_frame = tk.Frame(root)
        top_frame = tk.Frame(third_frame)
        # fourth_frame = tk.Frame(root)
        
        # Pack frames to control layout
        left_frame.pack(side=tk.LEFT, padx=10, pady=10)
        right_frame.pack(side=tk.LEFT, padx=10, pady=10)
        third_frame.pack(side=tk.LEFT,padx=10, pady=10)
        top_frame.pack(side=tk.TOP,fill=tk.X,padx=10)
        # fourth_frame.pack(side=tk.LEFT,padx=10, pady=10)
        
        # Scan Index Controls in right_frame
        self.scan_label = tk.Label(right_frame, text=f"Scan: {self.scan}", font=bold_font)
        self.scan_label.pack(pady=5)  # Adds a little space above and below
        
        self.scan_button_inc = tk.Button(right_frame, text="Increase Scan Index", command=self.inc_scan)
        self.scan_button_inc.pack()
        
        self.scan_button_dec = tk.Button(right_frame, text="Decrease Scan Index", command=self.dec_scan)
        self.scan_button_dec.pack()
        
        # Frame Controls in left_frame
        self.frm_label = tk.Label(left_frame, text=f"Frame: {self.frm}", font=bold_font)
        self.frm_label.pack(pady=5)  # Adds a little space above and below
        
        self.frame_button_inc = tk.Button(left_frame, text="Increase Frame", command=self.inc_frame)
        self.frame_button_inc.pack()
        
        self.frame_button_dec = tk.Button(left_frame, text="Decrease Frame", command=self.dec_frame)
        self.frame_button_dec.pack()
        
        # Save Truth Button
        self.save_button = tk.Button(top_frame, text="Save Truth", command=self.save_truth)
        self.save_button.pack(side=tk.LEFT, padx=10, pady=10)
        
        # Clear Truth Button
        self.clear_button = tk.Button(top_frame, text="Clear Truth", command=self.clear_truth)
        self.clear_button.pack(side=tk.LEFT, padx=10, pady=10)
        
        # Box Coordinates Label at the bottom
        self.box_coords_label = tk.Label(third_frame, text="Box Coordinates: (eta1, tth1) - (eta2, tth2)")
        self.box_coords_label.pack(side=tk.TOP)
        
        self.preload_button = tk.Button(root, text="Preload ROIs", command=self.preload_ROIs)
        self.preload_button.pack(side=tk.LEFT, padx=10, pady=10)
        
        self.rect_selector = RectangleSelector(self.ax, self.onselect,
                                               useblit=True,
                                               button=[1],  # Left mouse button
                                               minspanx=5, minspany=5,
                                               spancoords='pixels', interactive=True)

        self.rect_selector.set_active(True)
        
        # self.scan_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        # self.frm_label.grid(row=0, column=1, padx=10, pady=10, sticky="w")
        # self.scan_button_inc.grid(row=1, column=0, padx=10, pady=10, sticky="w")
        # self.scan_button_dec.grid(row=1, column=1, padx=10, pady=10, sticky="w")
        # self.save_button.grid(row=5, column=0, columnspan=2, pady=10)
        
        # Keybinds  
        self.root.bind('<Left>', self.decrement_scan_key)
        self.root.bind('<Right>', self.increment_scan_key)
        self.root.bind('<Up>', self.decrement_frame_key)
        self.root.bind('<Down>', self.increment_frame_key)
        self.root.bind('a', self.decrement_scan_key)
        self.root.bind('d', self.increment_scan_key)
        self.root.bind('w', self.decrement_frame_key)
        self.root.bind('s', self.increment_frame_key)
        self.root.bind('q', self.save_truth_key)
        self.root.bind('1', self.dec_etaRoi_key)
        self.root.bind('3', self.inc_etaRoi_key)
        self.root.bind('2', self.reset_etaRoi_key)
        
        self.updatePlot()
        self.ax.set_title(f'Spot {spotInd}', fontsize=16, fontweight="bold")
        self.root.mainloop() 
         
    def reset_etaRoi_key(self,event):
        self.reset_etaRoi()
    def inc_etaRoi_key(self,event):
        self.inc_etaRoi()
    def dec_etaRoi_key(self,event):
        self.dec_etaRoi()
    def increment_scan_key(self, event):
        self.inc_scan()
    def decrement_scan_key(self, event):
        self.dec_scan()
    def increment_frame_key(self, event):
        self.inc_frame()
    def decrement_frame_key(self, event):
        self.dec_frame()
    def save_truth_key(self, event):
        self.save_truth()
        
    def inc_etaRoi(self):
        roiSize = self.params['roiSize']
        detectDist, mmPerPixel, ff_trans, ff_tilt = dt.loadYamlData(self.params,self.tthRoi,self.etaRoi)
        rad_dom, eta_dom = dt.polarDomain(detectDist,mmPerPixel,self.tthRoi,self.etaRoi,roiSize)
        self.etaRoi = eta_dom[-1]
        self.updatePlot()
    def dec_etaRoi(self):
        roiSize = self.params['roiSize']
        detectDist, mmPerPixel, ff_trans, ff_tilt = dt.loadYamlData(self.params,self.tthRoi,self.etaRoi)
        rad_dom, eta_dom = dt.polarDomain(detectDist,mmPerPixel,self.tthRoi,self.etaRoi,roiSize)
        self.etaRoi = eta_dom[0]
        self.updatePlot()
    def reset_etaRoi(self):
        self.etaRoi = self.eta0
        self.updatePlot()
    def inc_scan(self):  
        self.scan_ind = self.scan_ind + 1
        if self.scan_ind >= len(self.scanRange):
            self.scan_ind = self.scan_ind - 1
        elif self.scan_ind < 0:
            self.scan_ind = 0;
        self.scan = self.scanRange[self.scan_ind]
        self.updatePlot()
    def dec_scan(self):
        self.scan_ind = self.scan_ind - 1
        if self.scan_ind >= len(self.scanRange):
            self.scan_ind = self.scan_ind - 1
        elif self.scan_ind < 0:
            self.scan_ind = 0;
        self.scan = self.scanRange[self.scan_ind]
        self.updatePlot()
    def inc_frame(self):
        self.frm = wrapFrame(self.frm + 1)
        self.updatePlot()
    def dec_frame(self):
        self.frm = wrapFrame(self.frm - 1)
        self.updatePlot()
    
    def updatePlot(self):
        self.ax.cla()
        
        # 1 Show ROI
        showROI(self.ax,self.dataPath, self.scan, self.frm,\
                self.tthRoi, self.etaRoi, self.params)
        showInitial(self.ax,self.etaRoi,self.tthRoi,self.eta0,self.tth0,self.params)
        # 2 Show the track and truth if there is one
        [self.trackFound, self.om_ind1] = checkTrack(self.track_data,self.scan_ind,\
                                           self.scan,self.frm)
        [self.truthFound, self.om_ind2] = checkTruth(self.truth_data,self.scan_ind,\
                                           self.scan,self.frm)      
        if self.trackFound:
            track = self.track_data[self.scan_ind][self.om_ind1].copy()
            showTrack(self.ax,track,self.etaRoi,self.tthRoi,\
                      self.eta0,self.tth0,self.params)
        if self.truthFound:
            print(f'Truth found: scan={self.scan},frm={self.frm}')
            self.truth = self.truth_data[self.scan_ind][self.om_ind2].copy()
            showTruth(self.ax,self.truth,self.etaRoi,self.tthRoi,self.params)
            self.select_truth()
        
        if hasattr(self,'truth'):
            self.truth['scan'] = self.scan
            self.truth['frm'] = self.frm
            self.select_truth()
        
        # 4 Produce plot
        self.scan_label.config(text=f"Scan: {self.scan}")
        self.scan_label.pack()
        self.frm_label.config(text=f"Frame: {self.frm}")
        self.frm_label.pack()
        plt.draw()
        self.canvas.draw()

    def onselect(self, eclick, erelease):
        # eclick and erelease are the points where the user clicked the mouse to start and end the selection
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        # Convert box to eta/tth values
        eta1,tth1,deta,dtth = pixToEtaTth(x1,y1,self.tthRoi,self.etaRoi,self.params)
        eta2,tth2,deta,dtth = pixToEtaTth(x2,y2,self.tthRoi,self.etaRoi,self.params)
        self.box_coords_label.config(text=f"Box: ({eta1:.4f}, {tth1:.4f}) - ({eta2:.4f}, {tth2:.4f})")
        self.truth = initTruth(eta1,eta2,tth1,tth2,self.scan,self.frm)

    def select_truth(self):
        eta1 = self.truth['eta1']
        eta2 = self.truth['eta2']
        tth1 = self.truth['tth1']
        tth2 = self.truth['tth2']
        x1,y1 = etaTthToPix(eta1,tth1,self.etaRoi,self.tthRoi,self.params)
        x2,y2 = etaTthToPix(eta2,tth2,self.etaRoi,self.tthRoi,self.params)
        self.box_coords_label.config(text=f"Box: ({eta1:.4f}, {tth1:.4f}) - ({eta2:.4f}, {tth2:.4f})")
        self.rect_selector.extents = (y1, y2, x1, x2)

    def save_truth(self):
        self.truth_data = addTruth(self.truth_data,self.scan_ind,self.truth)
        self.updatePlot()
        with open(self.truth_file, 'wb') as f:
            pickle.dump(self.truth_data, f)
            print('Saved truth data')
            
    def clear_truth(self):
        if self.truthFound:
            self.truth_data[self.scan_ind].pop(self.om_ind2)
            self.truthFound = False
        else:
            print('No truth to clear')
            
        self.updatePlot()
        with open(self.truth_file, 'wb') as f:
            pickle.dump(self.truth_data, f)
            print('Saved truth data')
        
    def preload_ROIs(self):
        print('Preloading Eiger data')
        for scan in self.scanRange:
            print(f'Loading scan {scan}')
            loadROI(self.dataPath,scan,self.frm,self.etaRoi,self.tthRoi,self.params)

def enterEtaTth(eta,tth):
    etaStr = input('Enter true Eta: ')
    try:
        newEta = float(etaStr)
    except:
        print('Keeping old eta')
        newEta = eta
    tthStr = input('Enter true Tth: ')
    try:
        newTth = float(tthStr)
    except:
        print('Keeping old tth')
        newTth = tth
    return newEta, newTth

def addTruth(truth_data,scan_ind,truth):
    # Make sure tracks are sort by omega when adding to list
    if len(truth_data[scan_ind]) == 0:
        truth_data[scan_ind].append(truth.copy())
    elif truth_data[scan_ind][0]['frm'] == truth['frm']:
        truth_data[scan_ind][0] = truth.copy()
    elif truth_data[scan_ind][0]['frm'] > truth['frm']:
        truth_data[scan_ind].insert(0,truth.copy())
    elif truth_data[scan_ind][-1]['frm'] < truth['frm']:
        truth_data[scan_ind].append(truth.copy())
    else:
        for i in range(len(truth_data[scan_ind])-1):
            if (truth_data[scan_ind][i]['frm'] < truth['frm']) and\
               (truth_data[scan_ind][i+1]['frm'] > truth['frm']):
                truth_data[scan_ind].insert(i+1,truth.copy())
                'Actually added?'
                break 
    return truth_data
    
def initTruth(eta1,eta2,tth1,tth2,scan,frm):      
    truth = {
            'eta1': eta1,
            'eta2': eta2,
            'tth1': tth1,
            'tth2':tth2,
            'scan': scan,
            'frm': frm
             }
    return truth