# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 22:18:31 2024

@author: dpqb1
"""

from hexrd import transforms
import os
import spotfetch as sf
import numpy as np

topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
yamlFile = os.path.join(topPath,'c103_eiger_calibration.yml')

detectDist, mmPerPixel, trans, tilt = sf.loadYamlDataEiger(yamlFile)

# Ryaw = np.array([[np.cos(tilt[0]),0,np.sin(tilt[0])],
#                  [0,1,0],
#                  [-np.sin(tilt[0]),0,np.cos(tilt[0])]])
# Rpit = np.array([[1,0,0],
#                  [0,np.cos(tilt[1]),-np.sin(tilt[1])],
#                  [0,np.sin(tilt[0]),np.cos(tilt[0])]])
# Rrol = np.array([[np.cos(tilt[2]),-np.sin(tilt[1]),0],
#                  [np.sin(tilt[2]),np.cos(tilt[2]),0],
#                  [0,0,1]])
# R1 = np.matmul(Ryaw,Rpit)

# R = np.matmul(R1,Rrol)

# rMat_d = Rh
# rMat_s = np.eye(3)
# rMat_c = np.eye(3)
# tVec_d = np.array(ff_trans)
# tVec_s = np.ones((3,1))
# tVec_c = np.ones((3,1))

angs = [0.2,0.1]
r = detectDist*np.tan(angs[1])/mmPerPixel
x = r*np.cos(angs[0])
y = -r*np.sin(angs[0])

out = sf.applyTransTilt(x,y,tilt,trans)

