# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 13:35:59 2017

@author: dbanco02
"""

## Init
from __future__ import division
from __future__ import print_function

import numpy as np
import os
import RingModel as RM
import EllipticModels as EM
import random
from functools import partial
import multiprocessing


ringModels_synth = np.load('ringModels_synth.npy')
print(len(ringModels_synth))
for i in range(len(ringModels_synth)):
    print(ringModels_synth[i].times)