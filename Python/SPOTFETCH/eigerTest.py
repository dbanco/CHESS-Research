# -*- coding: utf-8 -*-
"""
Created on Tue May 28 19:16:02 2024

@author: dpqb1
"""

import spotfetch as sf
import matplotlib.pyplot as plt


fname = "C:\\Users\\dpqb1\\Documents\\Data\\mruby-0401-1_0001_EIG16M_CdTe_000_008682_data_000001.h5"
# yamlFile = "C:\\Users\\dpqb1\\Documents\\Data\\mruby_eiger_calibration_single_grain_v01.yml"
# t = 1000
# t2 = 50
# img  = sf.load_eiger(fname, yamlFile, t)
# plt.imshow(img)

from hexrd import imageseries

ims = imageseries.open(fname, format='eiger-stream-v1')
imageseries.write(ims, 'output.npz', 'frame-cache', style='npz', threshold=200)