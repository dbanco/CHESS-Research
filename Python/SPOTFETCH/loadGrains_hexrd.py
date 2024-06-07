# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 10:28:22 2024

@author: dpqb1
"""
import yaml
from hexrd import instrument
from hexrd import imageseries
import numpy as np
import matplotlib.pyplot as plt
import os 

chunked_yaml_filepath = '/nfs/chess/aux/reduced_data/cycles/2023-2/id3a/shanks-3731-a/ti-2-dex-refinement-1.yml'
# ff1_filepath = '/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/1/ff/ff1_000106.h5'
# ff2_filepath = '/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/1/ff/ff2_000106.h5'

chunked_yaml_filepath = 'D:\\CHESS_raw_data\\ti-2-exsitu\\ti-2-dex-refinement-1.yml'
ff1_filepath = 'D:\\CHESS_raw_data\\ti-2-exsitu\\12\\ff\\ff1_000098.h5'
ff2_filepath = 'D:\\CHESS_raw_data\\ti-2-exsitu\\12\\ff\\ff2_000098.h5'

with open(chunked_yaml_filepath) as f:
    icfg = yaml.load(f, Loader=yaml.SafeLoader)
instr = instrument.HEDMInstrument(instrument_config=icfg)

panel_rois = {(det_key, panel.roi) for det_key, panel in instr.detectors.items()}

chunked_ims_dict = {}
for subpanel_key, roi in panel_rois:
    if subpanel_key.startswith('ff2'):
       ims = imageseries.open(ff2_filepath,format='hdf5',path='/imageseries')
       oplist = [('flip', 'h'),('rectangle', roi)]
    else:
        ims = imageseries.open(ff1_filepath,format='hdf5',path='/imageseries')
        oplist = [('flip', 'v'),('rectangle', roi)]
    chunked_ims_dict[subpanel_key] = imageseries.process.ProcessedImageSeries(ims, oplist)

img1 = chunked_ims_dict[subpanel_key][0]
plt.imshow(img1)
plt.clim(150,350)

full_roi = ((0,3072),(0,3888))
ims2 = imageseries.open(ff2_filepath,format='hdf5',path='/imageseries')
oplist = [('flip', 'h'),('rectangle',full_roi)]
r_ims = imageseries.process.ProcessedImageSeries(ims2, oplist)

ims1 = imageseries.open(ff1_filepath,format='hdf5',path='/imageseries')
oplist = [('flip', 'v'),('rectangle', full_roi)]
l_ims = imageseries.process.ProcessedImageSeries(ims1, oplist)
    
left = np.max(l_ims[0:3], axis=0)
right = np.max(r_ims[0:3], axis=0)
full_ims = np.concatenate((left,np.fliplr(right)),1)
b_median = np.median(full_ims)
full_ims = full_ims - 1.1 * b_median
full_ims[full_ims < 0] = 0

plt.imshow(np.flipud(full_ims))
plt.clim(0,350)

