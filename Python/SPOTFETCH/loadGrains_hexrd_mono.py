# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 11:21:02 2024

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
# chunked_yaml_filepath = 'D:\\CHESS_raw_data\\ti-2-exsitu\\ti-2-dex-refinement-1.yml'

monolithic_yaml_filepath = 'D:\\CHESS_raw_data\\ti-2-exsitu\\dex-refined-1.yml'
ff1_filepath = 'D:\\CHESS_raw_data\\ti-2-exsitu\\12\\ff\\ff1_000098.h5'
ff2_filepath = 'D:\\CHESS_raw_data\\ti-2-exsitu\\12\\ff\\ff2_000098.h5'


with open(monolithic_yaml_filepath) as f:
    icfg = yaml.load(f, Loader=yaml.SafeLoader)
instr = instrument.HEDMInstrument(instrument_config=icfg)

panel_rois = {(det_key, panel.roi) for det_key, panel in instr.detectors.items()}

monolithic_ims_dict = {}
for subpanel_key, roi in panel_rois:
    if subpanel_key == 'ff2':
        ims = imageseries.open(ff2_filepath,format='hdf5',path='/imageseries')
        oplist = [('flip', 'h')]
    else:
        ims = imageseries.open(ff1_filepath,format='hdf5',path='/imageseries')
        oplist = [('flip', 'v')]
    monolithic_ims_dict[subpanel_key] = imageseries.process.ProcessedImageSeries(ims, oplist)

# %%  
left = np.max(monolithic_ims_dict['ff2'][1345:1346], axis=0)
right = np.max(monolithic_ims_dict['ff1'][1345:1346], axis=0)
left = np.max(monolithic_ims_dict['ff2'][1345:1346], axis=0)
right = np.max(monolithic_ims_dict['ff1'][1345:1346], axis=0)
full_ims = np.concatenate((left,np.fliplr(right)),1)

plt.figure()
plt.imshow((full_ims))
plt.clim(313,500)

# From Sven's message
# [1345,1780:1830,940:1020]
# np.array(monolithic_ims_dict[subpanel_key])
# monolithic_ims_dict[subpanel_key][1345]



