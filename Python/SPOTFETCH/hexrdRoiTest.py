# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 22:11:18 2024

@author: dpqb1
"""
import numpy as np 
import os
from matplotlib import pyplot as plt
from hexrd import instrument
from hexrd import xrdutil
from hexrd import transforms
import glob

import spotfetch as sf

# LOCAL
topPath = r"C:\Users\dpqb1\Documents\Data\VD_sim_processing"
exsituPath = os.path.join(topPath,r"state_0\simulation\outputs\c103_polycrystal_sample_2_state_0_layer_1_output_data.npz")
dataFile = os.path.join(topPath,r"state_*\simulation\outputs\c103_polycrystal_sample_2_state_*_layer_1_output_data.npz")

# %% 2. Load in or collect spots data
# Spots for each state file
state = 0
# for state in range(5):
spotsDir = os.path.join(topPath,f'state_{state}','simulation','outputs')
spotsOut = os.path.join(topPath,f'state_{state}')
        
    # sf.collectSpotsData(spotsOut, spotsDir)
#
spotData = np.load(os.path.join(spotsOut,'spots.npz'))

# Detector and tracking parameters
params = {}
params['detector'] = 'eiger_sim'
params['peak_func'] = "gaussian"
params['imSize'] = (5000,5000)
params['yamlFile'] = os.path.join(topPath,'c103_eiger_calibration.yml') #mruby_0401_eiger_calibration
params['roiSize'] = [40,40]
params['gamma'] = [4,5,9,6] #[eta,tth,fwhm_eta,fwhm_tth]
params['pool'] = 16
params['parallelFlag'] = False

fullscanRange = np.arange(5)
trackPath = os.path.join(topPath,'outputs_test')
fname = os.path.join(topPath,r"state_0\simulation\outputs\c103_polycrystal_sample_2_state_0_layer_1_output_data.npz")

# configuration_filepath = r'C:\Users\dpqb1\CHESS-Research\Python\SPOTFETCH\VD_configuration.yml'
# config = Configuration.open_file(configuration_filepath)[0]
# # Load the instrument
# instr = Utilities.load_instrument(config)
# output = dict.fromkeys(instr.detectors)
# for detector_id, panel in instr.detectors.items()

yamlFile = params['yamlFile']
detectDist, mmPerPixel, trans, tilt = sf.loadYamlDataEiger(yamlFile)

simData = np.load(fname)       
shp = simData['shape']
img = np.zeros((shp[0],shp[1]))

eigDetector = instrument.detector.Detector(rows=shp[0],cols=shp[1],
                                            pixel_size=(0.075,0.075),
                                            tvec = trans,tilt=tilt)
instr_cfg = sf.read_yaml(yamlFile)
instr_cfg['detector'] = instr_cfg['detectors']['eiger']

frame_i = 2

rowD = simData[f'{frame_i-2}_row']
colD = simData[f'{frame_i-2}_col']
datD = simData[f'{frame_i-2}_data']

for i in range(len(rowD)):
    img[rowD[i],colD[i]] = datD[i]

b = img

fig = plt.figure()
# b[b>100] = 100
plt.imshow(b)
plt.clim(np.median(b),np.median(b)+20)
plt.colorbar()

spotInds = sf.findSpots(spotData, frm = frame_i)

for k in spotInds:
    tth = spotData['tths_pred'][k]
    eta = spotData['etas_pred'][k]
    ome = spotData['omes_pred'][k]
    
    allAngs = np.array([[tth,eta,ome]])
    rMat_d = transforms.xf.makeDetectorRotMat(tilt)
    rMat_c = np.eye(3)#transforms.xfcapi.makeRotMatOfExpMap()
    chi = 0
    tVec_d = np.array(trans)
    tVec_c = np.zeros((3,1))
    tVec_s = np.zeros((3,1))
    distortion = None #np.eye(3)
    
#     # Try to run sven's code and see what the variables end up being
    det_xy, rmats_s, on_plane = xrdutil._project_on_detector_plane(allAngs,
                                                                rMat_d, rMat_c, chi,
                                                                tVec_d, tVec_c, tVec_s,
                                                                distortion)
    
    x_pred = spotData['Xs'][k]
    y_pred = spotData['Ys'][k]
    alt_xy = np.array([[x_pred,y_pred]])
    
    x_meas = spotData['Xm'][k]
    y_meas = spotData['Ym'][k]
    alt_xy2 = np.array([[x_meas,y_meas]])
    
    det_xy = alt_xy2
    # print(det_xy)
    # print(x_pred)
    # print(y_pred)
    # print('-------')
    
    det_pix = det_xy.ravel()/mmPerPixel
    frame_i = spotData['ome_idxs'][k]
    plt.plot(det_pix[0]+ (shp[1]-1)/2,
              -det_pix[1]+ (shp[0]-1)/2 ,'x',markersize=12)
    # break
             

r = detectDist*np.tan(tth)/mmPerPixel
ang_pixel_size = 1/r
patches = xrdutil.make_reflection_patches(instr_cfg,
                            np.array([[tth,eta]]),
                            np.array([[ang_pixel_size,ang_pixel_size]]), omega=None,
                            tth_tol=0.2, eta_tol=1.0,
                            rmat_c=np.eye(3), tvec_c=np.zeros((3, 1)),
                            npdiv=1, quiet=True)



# strip relevant objects out of current patch
# for patch in patches:
#     vtx_angs, vtx_xy, conn, areas, xy_eval, ijs = patch

#     # Remap ijs to pull just the detector pixels - one to one
#     x_max = np.max(ijs[0])
#     x_min = np.min(ijs[0])
#     y_max = np.max(ijs[1])
#     y_min = np.min(ijs[1])
    
#     x_remap = np.arange(x_min,x_max)
#     y_remap = np.arange(y_min,y_max)
#     xs,ys = np.meshgrid(x_remap,y_remap)
#     ijs_remap = [None]*2
#     ijs_remap[0] = xs
#     ijs_remap[1] = ys
    
#     prows, pcols = areas.shape
#     nrm_fac = areas/float(native_area)
#     nrm_fac = nrm_fac / np.min(nrm_fac)
    
#     # grab hkl info
#     hkl = hkls_p[i_pt, :]
#     hkl_id = hkl_ids[i_pt]
    
#     # edge arrays
#     tth_edges = vtx_angs[0][0, :]
#     delta_tth = tth_edges[1] - tth_edges[0]
    
#     eta_edges = vtx_angs[1][:, 0]
#     delta_eta = eta_edges[1] - eta_edges[0]
    
#     # need to reshape eval pts for interpolation
#     xy_eval = np.vstack([xy_eval[0].flatten(),xy_eval[1].flatten()]).T
    
#     # the evaluation omegas;
#     # expand about the central value using tol vector
#     ome_eval = np.degrees(ang_centers[i_pt, 2]) + ome_del
    
#     # ???: vectorize the omega_to_frame function to avoid loop?
#     frame_indices = [ome_imgser.omega_to_frame(ome)[0] for ome in ome_eval]
    
# patch_data = panel.interpolate_bilinear(xy_eval,ome_imgser[i_frame],
#                                         pad_with_nans=False).reshape(prows, pcols)  # * nrm_fac
    # outputs are:
    #     (tth_vtx, eta_vtx),
    #     (x_vtx, y_vtx),
    #     connectivity,
    #     subpixel_areas,
    #     (x_center, y_center),
    #     (i_row, j_col)

    


