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
    def __init__(self, root,read_path):
        self.root = root
        self.root.title("Live Data Plotter")
        self.read_path = read_path
        # Explicitly create a figure and pass it to subplots
        self.fig = Figure(figsize=(5, 8))
        self.gs = gridspec.GridSpec(5, 1, height_ratios=[1, 1, 1, 1, 1])
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=45)

        # Create subplots within the same figure
        self.axes = [self.fig.add_subplot(self.gs[i]) for i in range(5)]

        self.trackData = []
        self.i = 1

        # Dropdown menu for selecting plot
        self.dropdown_menus = []
        dropdown_options = ["FWHM_omega", "FWHM_eta", "FWHM_tth","Mean_omega","Mean_eta", "Mean_tth",]
        skip = 9
        for i in range(5):
            var = tk.StringVar()
            var.set(dropdown_options[i])
            dropdown = ttk.Combobox(self.root, textvariable=var, values=dropdown_options,width=10)
            dropdown.grid(row=skip*i, column=1,pady=1)
            dropdown.bind('<<ComboboxSelected>>', self.update_plot)
            self.dropdown_menus.append((dropdown, var))
        
        self.update_plots()
        
    def update_plot(self, event=None):
        selected_plot_types = [dropdown_var.get() for _, dropdown_var in self.dropdown_menus]
        for i, (ax, selected_plot_type, spot_number) in enumerate(zip(self.axes, selected_plot_types)):
            ax.clear()
            # ax.set_title(f'{selected_plot_type} of Spot {spot_number}')
                # Handle different plot types based on selected option
            if selected_plot_type == "FWHM_omega":
                T = len(self.trackData)
                FWHMomega = np.zeros((T,1))
                scan = np.zeros((T,1))
                for t in range(T):
                    if self.trackData[t] != []:
                        if self.trackData[t] != []:
                            FWHMomega[t] = sf.estFWHMomega(self.trackData[t])
                            scan[t] = self.trackData[t][0]['scan']
                    else:
                        FWHMomega[t] = None
                        scan[t] = None
                ax.plot(scan,FWHMomega,'-o')
                ax.set_ylabel(r"$FWHM_\eta$")  # Set y-axis label
                ax.set_xlabel("Scan #")  # Set x-axis label
                pass
            elif selected_plot_type == "FWHM_eta":
                T = len(self.trackData)
                FWHMeta = np.zeros((T,1))
                scan = np.zeros((T,1))
                for t in range(T):
                    if self.trackData[t] != []:
                        if self.trackData[t] != []:
                            FWHMeta[t] = self.trackData[t][0]['p'][3]
                            scan[t] = self.trackData[t][0]['scan']
                    else:
                        FWHMeta[t] = None
                        scan[t] = None
                ax.plot(scan,FWHMeta,'-o')
                ax.set_ylabel(r"$FWHM_\eta$")  # Set y-axis label
                ax.set_xlabel("Scan #")  # Set x-axis label
                pass
            elif selected_plot_type == "FWHM_tth":
                # Your plot logic for FWHM_tth
                T = len(self.trackData)
                FWHMtth = np.zeros((T,1))
                scan = np.zeros((T,1))
                for t in range(T):
                    if self.trackData[t] != []:
                        if self.trackData[t] != []:
                            FWHMtth[t] = self.trackData[t][0]['p'][4]
                            scan[t] = self.trackData[t][0]['scan']
                    else:
                        FWHMtth[t] = None
                        scan[t] = None
                ax.plot(scan,FWHMtth,'-o')
                ax.set_ylabel(r"$FWHM_{2\theta}$")  # Set y-axis label
                ax.set_xlabel("Scan #")  # Set x-axis label
                pass
            elif selected_plot_type == "Mean_omega":
                # Your plot logic for Mean_eta
                T = len(self.trackData)
                MEANomega = np.zeros((T,1))
                scan = np.zeros((T,1))
                for t in range(T):
                    if self.trackData[t] != []:
                        if self.trackData[t] != []:
                            MEANomega[t] = sf.estMEANomega(self.trackData[t])
                            scan[t] = self.trackData[t][0]['scan']
                    else:
                        MEANomega[t] = None
                        scan[t] = None
                ax.plot(scan,MEANomega,'-o')
                ax.set_ylabel(r"$\mu_\eta$")  # Set y-axis label
                ax.set_xlabel("Scan #")  # Set x-axis label
                pass
            elif selected_plot_type == "Mean_eta":
                # Your plot logic for Mean_eta
                T = len(self.trackData)
                MEANeta = np.zeros((T,1))
                scan = np.zeros((T,1))
                for t in range(T):
                    if self.trackData[t] != []:
                        if self.trackData[t] != []:
                            MEANeta[t] = self.trackData[t][0]['eta']
                            scan[t] = self.trackData[t][0]['scan']
                    else:
                        MEANeta[t] = None
                        scan[t] = None
                ax.plot(scan,MEANeta,'-o')
                ax.set_ylabel(r"$\mu_\eta$")  # Set y-axis label
                ax.set_xlabel("Scan #")  # Set x-axis label
                pass
            elif selected_plot_type == "Mean_tth":
                # Your plot logic for Mean_tth
                T = len(self.trackData)
                MEANtth = np.zeros((T,1))
                scan = np.zeros((T,1))
                for t in range(T):
                    if self.trackData[t] != []:
                        if self.trackData[t] != []:
                            MEANtth[t] = self.trackData[t][0]['tth']
                            scan[t] = self.trackData[t][0]['scan']
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
                spot_numbers = [spot_entry.get() for spot_entry in self.spot_entries]
                for i in range(5):
                    try:
                        k = spot_numbers[i]
                        with open(os.path.join(self.read_path,f'trackData_{k}.pkl'), 'rb') as f:
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
    read_path1 = '/mnt/scratch/dbanco/....TBD'
    spotInds1 = np.arange(20)
    start_gui(read_path1,spotInds1)

