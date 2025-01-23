# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 13:34:37 2025

Script to load in .pkl file containing a single 
roiTensor with shape [M_tth,M_eta,M_ome,T]

The ROI is extracted at the same [tth,eta] location at each omega and each scan.
In the data I have provided, M_ome = 7 and the spot of interest should at the 
very least appear centered in frame roiTensor[:,:,3,0] 

HEXD spot data, shared here has no spot tracking information at all. It is possible
for the spot to exit the ROI either in eta or omega. Spot 4650 for example, the
first spot, does not appear to move.


@author: Daniel Banco
"""
import pickle
import os 
import matplotlib.pyplot as plt
import numpy as np
import spotfetch as sf

from math import sqrt
from skimage import data
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage.color import rgb2gray

from scipy.signal import convolve
from skimage.feature import hessian_matrix
from scipy.ndimage import label
from sklearn.mixture import BayesianGaussianMixture


def init3d(shape):
    shape = (2*shape[0]+1,2*shape[1]+1,2*shape[2]+1)
    x = np.zeros(shape)
    y = np.zeros(shape)
    z = np.zeros(shape)
    for i in range(-shape[0],shape[0]):
        x[i,:,:] = i
    for i in range(-shape[1],shape[1]):
        y[:,i,:] = i
    for i in range(-shape[2],shape[2]):
        z[:,:,i] = i
    return x,y,z

def gauss_kernel(shape,sigma):
    x,y,z = init3d(shape)
    G = (2*np.pi*sigma**2)**(3/2)*np.exp(-(x**2+y**2+z**2)/(2*sigma**2))
    return G

def DoG1(f,sigma,dsigma,gamma=2):
    g = (gauss_kernel(f.shape,sigma + dsigma) - 
         gauss_kernel(f.shape,sigma))/(sigma*dsigma)
    return sigma**(gamma-1)*convolve(f,g,'same')                                  

# Collapse this function for plotting
def plotROI(roiTensor,num_figs):
    T = roiTensor.shape[3]
    M_ome = roiTensor.shape[2]
    
    num_cols = int(np.ceil(T/num_figs))
    scanRange = np.arange(29)
    
    fig_list = []
    axes_list = []
    
    # Create a figure and a set of subplots
    for i in range(num_figs):
        fig, axes = plt.subplots(M_ome, num_cols, figsize=(20, 15))
        
        # Remove x and y ticks for clarity
        for ax_row in axes:
            for ax in ax_row:
                ax.set_xticks([])
                ax.set_yticks([])
            
        fig_list.append(fig)
        axes_list.append(axes)
        
        i1 = 0 + i*num_cols
        i2 = num_cols-1 + i*num_cols
        j1 = scanRange[i1]
        if i2 >= T:
            j2 = scanRange[-1]
        else:
            j2 = scanRange[i2]
        # Add common labels
        fig_list[i].text(0.04, 0.5, r'$\omega$ frame', va='center', rotation='vertical', fontsize=24)
        fig_list[i].text(0.5, 0.04, 'Scan #', ha='center', fontsize=24)
        fig_list[i].text(0.5, 0.95, f'Spot {spotInd}, Scans {j1}-{j2}', ha='center', fontsize=32)      
        fig_list[i].subplots_adjust(wspace=0.05, hspace=0.01)
        
    for scan_ind in range(T):
        for om_ind in range(M_ome):
            i = int(np.floor(scan_ind/num_cols))
            j = np.mod(scan_ind,num_cols)
            ax = axes_list[i][om_ind,j]
            scan = scanRange[scan_ind]
            
            # 1. Show ROI
            roi = roiTensor[:,:,om_ind,scan_ind]
            ax.imshow(roi)
            ax.text(1, 2, f'max: {roi.max():.2f}', color='white', fontsize=12, weight='bold')
          
            # Label plots
            if om_ind == M_ome-1:
                ax.set_xlabel(f'{scan}')
        

# Load and plot ROI
spotInd = 4653
topDir = r'E:\Data\c103_processing\roiTensors_grain_44'
roiFile = os.path.join(topDir,f'roiTensor_{spotInd}.pkl')
with open(roiFile, 'rb') as f:
    roiTensor = pickle.load(f)
    

# Try 3d blob detection with DoG

roi = roiTensor[:,:,:,0]
x = roi.copy()
x[x<5] = 0

blobs, num_blobs, hess_mat = sf.detectBlobDoG(x)
RT_T,ST_T,AT_T = sf.blobFeaturesDoG(x,blobs,num_blobs,hess_mat)

X = np.transpose(np.array([RT_T,ST_T,AT_T]))

# Fit a Variational Bayesian Gaussian Mixture Model
vbgmm = BayesianGaussianMixture(
    n_components=3,  # Maximum number of components
    covariance_type='full',  # Full covariance matrices
    weight_concentration_prior_type='dirichlet_process',  # Non-parametric DP
    random_state=42
)

vbgmm.fit(X)

# Predict cluster labels
labels = vbgmm.predict(X)



# hessian = np.zeros((roi.shape[0]-2,roi.shape[1]-2,roi.shape[2]-2,3,3)) 
# for i in range(3):
#     for j in range(3):
#         h_temp = np.diff(np.diff(dog_norm,1,i),1,j)
#         hessian[:,:,:,i,j] = h_temp[:28,:28,:5]

# determin = np.zeros(hessian.shape[:3])
# for i1 in range(hessian.shape[0]):
#     for i2 in range(hessian.shape[1]):
#         for i3 in range(hessian.shape[2]):
#             determin[i1,i2,i3] = np.linalg.det(hessian[i1,i2,i3,:,:])


# Detect blobs across all omegas and scans                
for t in range(1):
    M_ome = 7
    blobs_list = []
    fig, axes = plt.subplots(3, 3, figsize=(9, 3), sharex=True, sharey=True)
    fig2, axes2 = plt.subplots(3, 3, figsize=(9, 3), sharex=True, sharey=True)
    ax = axes.ravel()
    ax2 = axes2.ravel()
    
    blob_1 = blobs.copy()
    blob_1[blobs != 1] = 0
    
    for j in range(M_ome):
        image = roiTensor[:,:,j,t]
        pimg = blobs[:,:,j]
        
        blobs_dog = blob_dog(image, max_sigma=8, threshold=5)
        blobs_dog[:, 2] = blobs_dog[:, 2] * sqrt(2)
        
        blobs_list.append(blobs_dog)
        
        ax[j].imshow(image)
        ax2[j].imshow(pimg)
        # for blob in blobs_list[j]:
        #     y, x, r = blob
        #     c = plt.Circle((x, y), r, color='red', linewidth=2, fill=False)
        #     ax[j].add_patch(c)
        ax[j].set_axis_off()
        
    plt.tight_layout()
    plt.show()

# from scipy.ndimage import label, find_objects
# roi_t0 = roiTensor[:,:,:,0]
# bin_img = roi_t0 >5
# l_img, n_labels= label(bin_img)

# max_label = 0
# max_count = 0
# for lab in range(1,n_labels):
#     count = np.sum(l_img == lab)
#     print(f'Label {lab}: {count}')
#     if count > max_count:
#         max_label = lab
#         max_count = count

# bin_mask = np.array(l_img == max_label)

# fig, axes = plt.subplots(3, 3, figsize=(9, 9), sharex=True, sharey=True)
# ax = axes.ravel()

# M_ome = 7
# for j in range(M_ome):
#     masked_img = np.multiply(roi_t0[:,:,j],bin_mask[:,:,j]+1-1)
#     ax[j].imshow(masked_img)
    
# objs = find_objects(l_img)

# plt.show(l_img)
# num_figs = 2 # Number figures across which to show roi data
# plotROI(roiTensor,num_figs)

