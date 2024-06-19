# -*- coding: utf-8 -*-
"""
Created on Tue May 28 19:16:02 2024

@author: dpqb1
"""

import spotfetch as sf
import matplotlib.pyplot as plt
import h5py
import hdf5plugin
from hexrd import imageseries
import numpy as np


fname = '/mnt/scratch/pagan-dwellfatigue-series-cycle-2024-2/ti6242-sample1-2_0011_EIG16M_CdTe_000001.h5'
yamlFile = '/mnt/scratch/dbanco/mruby_eiger_calibration_single_grain_v01.yml'
t = 10
t2 = 50

ims = imageseries.open(fname, format='eiger-stream-v1')
img = ims[0,:,:].copy()
mask = np.where(img==img.max())
img[mask] = 0
thresh = 50
img[img>50] = 50
plt.imshow(img)
plt.colorbar()




# file1 = h5py.File(fname, 'r')
# print(file1.keys())
# group = file1['data']
# print(group.keys())
# group2 = file1['metadata']
# meta1 = group2['output_metadata']
# data1 = group['2']
# data2 = data1['data']
# print(data2.shape)

# plt.imshow(img)

# img  = sf.load_eiger(fname, yamlFile, t)
# plt.imshow(img)

# from hexrd import imageseries

# ims = imageseries.open(fname, format='eiger-stream-v1')
# imageseries.write(ims, 'output.npz', 'frame-cache', style='npz', threshold=200)