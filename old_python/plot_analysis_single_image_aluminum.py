# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:01:01 2016

Batch peak fitting

@author: Dan
"""

## Init
import numpy as np
from numpy.linalg import norm 
import matplotlib.pyplot as plt
import matplotlib
import os
import time

from scipy.signal import argrelmax
from sklearn.linear_model import Lasso

import DataReader
import RingImageProcessing as RingIP
import EllipticModels as EM
import PeakFitting as peak
import DataAnalysis as DA
reload(RingIP)
reload(DataReader)
### Data directory
data_dir = os.path.join('F:','CHESS_raw_data')
out_dir  = os.path.join(data_dir,'out')

#rsf_out = np.load('results_exp_final_radial_sum_fixed.npy')
rsf_out = np.load('results_exp_aluminum_fixed.npy')

specimen_name         = 'al7075_mlf'
step_names            = ['initial',  '1turn',    '2turn',    '3turn',    'unload']
dic_files             = ['dic_4536', 'dic_4537', 'dic_4538', 'dic_4539', 'dic_4540']
dark_dirs             = [ 68,         277,        485,        692,        899]
init_dirs             = [ 68,         277,        485,        692,        899]
dic_center            = [0.16, 2.9]       # sample center (mm) in vic2d coordinates 
x_range, x_num        = [-6, 6], 5         # x-ray measurements, range in mm
y_range, y_num        = [-5, 5], 41        # x-ray measurements, range in mm
detector_dist         = 3289.95            # pixels
true_center           = [1020.67, 1024.61] # [row, column] of detector image center in pixels (shifted by 1 for python index)
e_rng                 = [-0.012, 0.012]    # elastic strain range
p_rng                 = [-0.024, 0.024]    # plastic strain range
t_rng                 = [-0.036, 0.036]    # total strain range
E, G, v               =  71.7, 26.9, 0.33  # elastic modulus (GPa), shear modulus (GPa), poisson's ratio

ring_name             = 'al_311'
radius                = 370                # ring radius in pixels
dr                    = 30                 # half of ring width in pixels
min_amp               = 25                 # minimum acceptable peak amplitude
vec_frac              = 0.25               # fraction of peaks that must be acceptable


sample                = DA.Specimen(specimen_name, data_dir, out_dir, 
                                    step_names, dic_files, dark_dirs, 
                                    init_dirs, dic_center, x_range, x_num, 
                                    y_range, y_num, detector_dist, true_center, 
                                    e_rng, p_rng, t_rng, E, G, v)      


num_theta = 2400
dtheta = 2*np.pi/num_theta
num_var = 100


#%% Image select   
img_num = 75
load_i = 0
rsf_ex = rsf_out[load_i][img_num]

#%% View single ring fit and image 
for i in range(5):
    img = sample.load_image(img_num,i) 
    plt.figure(i)
    plt.imshow(img,cmap='jet',vmin=0,vmax=200, interpolation='nearest')

a,b,xc,zc = RingIP.ring_fit_nonlin(img,370,30,pf1=1)


#%% View multiple ellipse fits
l_steps = [0,1,2,3,4]
img_num = 75
for i, load_i in enumerate(l_steps):
        plt.figure(load_i)
        rsf_ex = rsf_out[load_i][img_num]
        img = sample.load_image(img_num,load_i)
        rsf_ex.plot_ellipse(img)  

#%% Plot ellipse fit and interpolation points
reload(RingIP)
# Load ring image
ring_image = sample.load_image(img_num,load_i) 

# Fit ellipse to ring
a,b,xc,zc = RingIP.ring_fit_nonlin(ring_image,370,20,true_center,4,pf1=True)

radial_distance = 8 
theta_domain = np.linspace(0,2*np.pi-dtheta,num_theta)
# Azimuthal projection with radial sum
f_az = RingIP.azimuthal_projection_ellipse_radial_sum(ring_image,
                                              [xc,zc],a,b,
                                              theta_domain,
                                              radial_distance,pf1=4)

check = f_az[0] - rsf_ex.f

plt.figure(11)
plt.plot(theta_domain,f_az[0],'bo-')
plt.figure(12)
plt.plot(theta_domain,rsf_ex.f,'bo-')

#%%


plt.figure(3)
plt.rc('text',usetex=True)
#plt.rc('font',family='serif')


font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)
plt.xlabel(r'Angle $\eta$ (radians)')
plt.ylabel(r'Intensity $f(\eta)$')
plt.tick_params(axis='both', which='major', labelsize=16)
plt.plot(rsf_ex.theta,rsf_ex.f,'b-')

plt.figure(4)
plt.imshow(img,cmap='jet',vmin=0,vmax=200, interpolation='nearest')
xx, yy = RingIP.gen_ellipse(rsf_ex.a,rsf_ex.b,rsf_ex.xc,rsf_ex.zc,0,rsf_ex.theta)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.plot(xx,yy,'wo-')

plt.figure(5)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.plot(rsf_ex.theta,rsf_ex.f,'ob',label='Data')
plt.plot(rsf_ex.theta,rsf_ex.fit,'-r',label='Fit')
plt.xlabel(r'Angle $\eta$ (radians)')
plt.ylabel(r'Intensity $f(\eta)$')
plt.legend(fontsize=14)

plt.figure(6)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.plot(rsf_ex.theta,rsf_ex.f,'ob',label='Data')
rsf_ex.plot_bases_wrap(dtheta,6)
plt.xlabel(r'Angle $\eta$ (radians)')
plt.ylabel(r'Intensity $f(\eta)$')

index_max = np.argmax(rsf_ex.coefs)

print(rsf_ex.theta[int(rsf_ex.means[index_max])])
len(argrelmax(rsf_ex.f)[0])

#%% Plot fit but try to remove DC component

plt.figure(1)        
plt.ylabel('Intensity')
plt.xlabel('Angle (Radians)')

plt.plot(rsf_ex.theta,rsf_ex.f,'bo-') 
 
# Construct B
dom = np.arange(len(rsf_ex.coefs))
keep = dom[np.abs(rsf_ex.coefs)>0]
full_fit = np.zeros(len(rsf_ex.theta))

for i in keep:
    if(rsf_ex.coefs[i]<10):
        basis = EM.gaussian_basis_wrap(len(rsf_ex.f),dtheta,
                                  rsf_ex.means[i],rsf_ex.variances[i])
        plt.figure(1)
        plt.plot(rsf_ex.theta,rsf_ex.coefs[i]*basis,'r-') 
    
        full_fit += rsf_ex.coefs[i]*basis
plt.plot(rsf_ex.theta,full_fit,'g-')


#%% View an individual histogram and fit plots
num_theta = 2400
dtheta = 2*np.pi/num_theta
num_bins = 30

img_num = 13
load_steps = [0,4]
bar_colors = ['r','b','g','m','k']

hist_init, bins_init = np.histogram(rsf_out[0][0].variances,bins=num_bins,normed=1)
bar_width = (bins_init[1]-bins_init[0])/(2*len(load_steps))


for i, load_step in enumerate(load_steps):
    rsf_ex = rsf_out[load_step][img_num]
    indx = rsf_ex.coefs > 0
    variances = [rsf_ex.variances[indx]]
    print(np.mean(variances))
    hist_ex,_ = np.histogram(variances,bins=bins_init,normed=1)

    plt.figure(0)
    plt.bar(bins_init[0:-1]+i*bar_width,hist_ex,width=bar_width,color=bar_colors[i],label=str(load_step))
    plt.legend()

    plt.figure(1)
    plt.subplot(len(load_steps),1,i+1)
    norm_const = np.trapz(rsf_ex.f,dx=dtheta)
    plt.plot(rsf_ex.f/norm_const,'b-')
    
#%% Let's read in the strain data
components = ['exx','eyy','exy']
load_steps = ['initial','1turn','2turn','3turn','unload']    
    
base_path = 'F:CHESS\\out\\al7075_mlf\\al_311\\strain\\al_311_plastic_' 
import csv

strain = np.zeros((3,5,205))
for i_comp, comp in enumerate(components):
    for i_load, load_name in enumerate(load_steps):
        path = base_path + comp +'_' + load_name + '.txt'
        
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
            for i, row in enumerate(reader):
                strain[i_comp,i_load,i] =  float(row[2])
    

