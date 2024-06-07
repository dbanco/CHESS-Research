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
import numpy as np
import pickle
import threading
import time
import os

class DataPlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("Live Data Plotter")
        self.read_path = "C:\\Users\\dpqb1\\Documents\\Data\\indexed_grains\\outputs"
        self.cb_present = [False, False, False, False, False]
        self.cb = [None, None, None, None, None]
        # Explicitly create a figure and pass it to subplots
        self.fig = Figure(figsize=(5, 10))
        self.gs = gridspec.GridSpec(5, 1, height_ratios=[1, 1, 1, 1, 1])
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=45)

        # Create subplots within the same figure
        self.axes = [self.fig.add_subplot(self.gs[i]) for i in range(5)]

        self.trackData = []
        self.i = 1

        # Dropdown menu for selecting plot
        self.dropdown_menus = []
        dropdown_options = ["ROI/Track", "Omega Detections", "FWHM_eta", "FWHM_tth", "Mean_eta", "Mean_tth"]
        skip = 9
        for i in range(5):
            var = tk.StringVar()
            var.set(dropdown_options[i])
            dropdown = ttk.Combobox(self.root, textvariable=var, values=dropdown_options,width=10)
            dropdown.grid(row=skip*i, column=1,pady=1)
            dropdown.bind('<<ComboboxSelected>>', self.update_plot)
            self.dropdown_menus.append((dropdown, var))
        
        # Entry fields for entering spot numbers
        self.spot_entries = []
        for i in range(5):
            spot_entry = ttk.Entry(self.root,width=3)
            spot_entry_label = ttk.Label(self.root, text="Spot #:")
            spot_entry_label.grid(row=skip*i+1, column=1,sticky="W",padx=5)
            spot_entry.grid(row=skip*i+1, column=1,padx=5,pady=1)
            spot_entry.insert(0,"0")
            self.spot_entries.append(spot_entry)
 
        # Entry field for updating all spot numbers
        self.update_all_label = ttk.Label(self.root, text="Change All Spot #s:")
        self.update_all_label.grid(row=44, column=1, sticky="W", padx=5)
        self.update_all_entry = ttk.Entry(self.root, width=4)
        self.update_all_entry.grid(row=45, column=1, padx=5, pady=5)
        self.update_all_entry.bind("<Return>", self.update_all_spot_numbers)

    
        # Button to trigger plotting
        self.plot_button = ttk.Button(self.root, text="Plot", command=self.plot)
        self.plot_button.grid(row=45, column=0,columnspan=2,pady=10)
        
        self.update_plots()
        
    def update_all_spot_numbers(self, event):
        new_spot_number = self.update_all_entry.get()
        for spot_entry in self.spot_entries:
            spot_entry.delete(0, tk.END)
            spot_entry.insert(0, new_spot_number)
            
    def update_plot(self, event=None):
        selected_plot_types = [dropdown_var.get() for _, dropdown_var in self.dropdown_menus]
        spot_numbers = [spot_entry.get() for spot_entry in self.spot_entries]
        for i, (ax, selected_plot_type, spot_number) in enumerate(zip(self.axes, selected_plot_types, spot_numbers)):
            ax.clear()
            # ax.set_title(f'{selected_plot_type} of Spot {spot_number}')
            
            
            if spot_number.isnumeric():
                k = int(spot_number)
            
                # Handle different plot types based on selected option
                if selected_plot_type == "ROI/Track":
                    # Your plot logic for ROI/Track
                    T = len(self.trackData)
                    numROIs = 4
                    full_roi = self.trackData[T-numROIs-1][k][0]['roi']
                    maxVal = max(full_roi.ravel())
                    for t in range(T-numROIs,T-1):
                    # if self.trackData[T-2] != []:
                        full_roi = np.hstack((full_roi,maxVal*np.ones((40,2)),self.trackData[T-2][k][0]['roi']))
                    
                    img=ax.imshow(full_roi)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    if not self.cb[i] == None:
                        self.cb[i].update_normal(img)
                    else:
                        self.cb[i] = self.fig.colorbar(img, ax=ax)

                    # Calculate the starting and ending indices for rows and columns
                    boxSize = 10
                    tt = 0
                    for t in range(T-numROIs-1,T-1):
                        x = self.trackData[t][k][0]['p'][1]
                        y = self.trackData[t][k][0]['p'][2]
                        start_col = round(x - boxSize//2) + 42*tt
                        start_row = round(y - boxSize//2)
                        tt += 1
                        rect = plt.Rectangle((start_col, start_row), boxSize, boxSize,
                                                  linewidth=1, edgecolor='r', facecolor='none')
                        ax.add_patch(rect)
                        
                    pass
                elif selected_plot_type == "Omega Detections":
                    # Your plot logic for Omega Detections
                    T = len(self.trackData)
                    for t in range(T):
                        if self.trackData[t] != []:
                            L = len(self.trackData[t][k])
                            if L > 0:
                                for j in range(L):
                                    omega = sf.frameToOmega(self.trackData[t][k][j]['frm'])
                                    scan = self.trackData[t][k][j]['scan']
                                    ax.scatter(scan,omega,marker='s',color='b')
                                    ax.set_ylabel(r"$\omega$")  # Set y-axis label
                                    ax.set_xlabel("Scan #")  # Set x-axis label
                            
                    pass
                elif selected_plot_type == "FWHM_eta":
                    # Your plot logic for FWHM_eta
                    T = len(self.trackData)
                    FWHMeta = np.zeros((T,1))
                    scan = np.zeros((T,1))
                    for t in range(T):
                        if self.trackData[t] != []:
                            if self.trackData[t][k] != []:
                                FWHMeta[t] = self.trackData[t][k][0]['p'][3]
                                scan[t] = self.trackData[t][k][0]['scan']
                        else:
                            FWHMeta[t] = None
                            scan[t] = None
                    ax.plot(scan,FWHMeta,'-o')
                    ax.set_ylabel(r"$FWHM_\eta$")  # Set y-axis label
                    ax.set_xlabel("Scan #")  # Set x-axis label
                    pass
                elif selected_plot_type == "Mean_eta":
                    # Your plot logic for Mean_eta
                    T = len(self.trackData)
                    MEANeta = np.zeros((T,1))
                    scan = np.zeros((T,1))
                    for t in range(T):
                        if self.trackData[t] != []:
                            if self.trackData[t][k] != []:
                                MEANeta[t] = self.trackData[t][k][0]['eta']
                                scan[t] = self.trackData[t][k][0]['scan']
                        else:
                            MEANeta[t] = None
                            scan[t] = None
                    ax.plot(scan,MEANeta,'-o')
                    ax.set_ylabel(r"$\mu_\eta$")  # Set y-axis label
                    ax.set_xlabel("Scan #")  # Set x-axis label
                    pass
                elif selected_plot_type == "FWHM_tth":
                    # Your plot logic for FWHM_tth
                    T = len(self.trackData)
                    FWHMtth = np.zeros((T,1))
                    scan = np.zeros((T,1))
                    for t in range(T):
                        if self.trackData[t] != []:
                            if self.trackData[t][k] != []:
                                FWHMtth[t] = self.trackData[t][k][0]['p'][4]
                                scan[t] = self.trackData[t][k][0]['scan']
                        else:
                            FWHMtth[t] = None
                            scan[t] = None
                    ax.plot(scan,FWHMtth,'-o')
                    ax.set_ylabel(r"$FWHM_{2\theta}$")  # Set y-axis label
                    ax.set_xlabel("Scan #")  # Set x-axis label
                    pass
                elif selected_plot_type == "Mean_tth":
                    # Your plot logic for Mean_tth
                    T = len(self.trackData)
                    MEANtth = np.zeros((T,1))
                    scan = np.zeros((T,1))
                    for t in range(T):
                        if self.trackData[t] != []:
                            if self.trackData[t][k] != []:
                                MEANtth[t] = self.trackData[t][k][0]['tth']
                                scan[t] = self.trackData[t][k][0]['scan']
                        else:
                            MEANtth[t] = None
                            scan[t] = None
                    ax.plot(scan,MEANtth,'-o')
                    ax.set_ylabel(r"$\mu_{2\theta}$")  # Set y-axis label
                    ax.set_xlabel("Scan #")  # Set x-axis label
                    pass
                
        self.fig.tight_layout(pad=3.0)
        self.canvas.draw()
            
    def plot(self):
        self.update_plot()    
    
    def update_plots(self):
        def read_data():
            while True:
                try:
                    with open(os.path.join(self.read_path,'trackData.pkl'), 'rb') as f:
                        self.trackData = pickle.load(f)
                    self.root.after(0, self.update_plot)  # Schedule update_plot to run in the main thread
                except FileNotFoundError:
                    self.trackData = []
                except pickle.UnpicklingError:
                    self.trackData = []
                time.sleep(1)
        
        threading.Thread(target=read_data, daemon=True).start()
        
def start_gui():
    root = tk.Tk()
    app = DataPlotter(root)
    root.mainloop()

if __name__ == "__main__":
    start_gui()

