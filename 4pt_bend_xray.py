"""
analysis for Chess beam run

@author: Kenny Swartz
07/12/2016
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import os, importlib
os.chdir('/home/kswartz92/Documents/python_scripts/local_chess_files')
import DataReader   as DR
import DataAnalysis as DA
DR       = importlib.reload(DR)
DA       = importlib.reload(DA)
data_dir = '/media/kswartz92/Swartz/chess/'
out_dir  = '/home/kswartz92/Documents/chess_xray_data/'
#%%
specimen_name         = 'ti64_notched'
x_range               = [-6,6]          # mm
x_num                 = 13
z_range               = [-4,4]          # mm
z_num                 = 41
dic_center            = [-0.42, 2.236]  # sample center in vic2d coordinates(mm)
step_names            = ['initial',  '1000N',    '2000N',    'unload']
dic_files             = ['dic_2475', 'dic_2476', 'dic_2477', 'dic_2478']
dark_dirs             = [ 5,          544,        1081,       1618]
init_dirs             = [ 6,          545,        1082,       1619]
detector_dist         = 4261.46         # pixels
true_center           = [1026, 1028]    # [row, column] of center in pixels
e_rng                 = [-0.012, 0.012]
p_rng                 = [-0.036, 0.036]
t_rng                 = [-0.048, 0.048]
#%%    alpha phase
ring_name             = 'ti_100'
left                  = [508, 548]
right                 = [1506, 1546]
top                   = [503, 543]
bottom                = [1501, 1541]
lambda_h              = 0.4
lambda_v              = 0.1 
#%%    beta phase
ring_name             = 'ti_200'
left                  = [215, 245]
right                 = [1810, 1840]
top                   = [211, 241]
bottom                = [1804, 1834]
lambda_h              = 0.2
lambda_v              = 0.2
#%%
specimen_name         = 'al7075_mlf'
x_range               = [-6,6]          # mm
x_num                 = 5
z_range               = [-5,5]          # mm
z_num                 = 41
dic_center            = [0.16, 2.75]    # sample center in vic2d coordinates(mm)
step_names            = ['initial',  '1turn',    '2turn',    '3turn',    'unload']
dic_files             = ['dic_4536', 'dic_4537', 'dic_4538', 'dic_4539', 'dic_4540']
dark_dirs             = [ 68,         277,        485,        692,        899]
init_dirs             = [ 68,         277,        485,        692,        899]
detector_dist         = 3293.09         # pixels
true_center           = [1027, 1027]    # [row, column] of center in pixels
e_rng                 = [-0.012, 0.012]
p_rng                 = [-0.024, 0.024]
t_rng                 = [-0.036, 0.036]

ring_name             = 'al_311'
left                  = [278, 338]
right                 = [1713, 1773]
top                   = [274, 334]
bottom                = [1709, 1769]
lambda_h              = 2.0
lambda_v              = 0.8
#%%
#    create sample and ring objects
sample                = DA.Specimen(specimen_name, data_dir, out_dir, x_range, x_num, z_range, z_num, step_names, dark_dirs, init_dirs, detector_dist, true_center)
ring                  = DA.Ring(ring_name, sample, left, right, top, bottom, lambda_h, lambda_v)

# create appropriate coordinate arrays

if sample.name == 'ti64_notched' or specimen_name=='ti64_plain':
    xa                    = np.linspace(sample.x_range[1], sample.x_range[0], num=x_num, endpoint=True)
    za                    = np.linspace(sample.z_range[1], sample.z_range[0], num=z_num, endpoint=True)
    z2d, x2d              = np.meshgrid(za, xa)
if sample.name == 'al7075_plain':
    xa                    = np.linspace(x_range[1], x_range[0], num=x_num, endpoint=True)
    za                    = np.linspace(z_range[1], z_range[0], num=z_num, endpoint=True)
    x2d, z2d              = np.meshgrid(xa, za)
if sample.name == 'al7075_mlf':
    xa                    = np.linspace(x_range[0], x_range[1], num=x_num, endpoint=True)
    za                    = np.linspace(z_range[0], z_range[1], num=z_num, endpoint=True)
    x2d, z2d              = np.meshgrid(xa, za)

x1d, z1d                  = x2d.flatten(), z2d.flatten()

#    read peak diameters if they have been fit, if not fit here

orient                = 'h'
try:
    l_cent, l_errs, l_amps   = np.zeros((sample.n_load_step, sample.n_data_pt)), np.zeros((sample.n_load_step, sample.n_data_pt)), np.zeros((sample.n_load_step, sample.n_data_pt))
    u_cent, u_errs, u_amps   = np.zeros((sample.n_load_step, sample.n_data_pt)), np.zeros((sample.n_load_step, sample.n_data_pt)), np.zeros((sample.n_load_step, sample.n_data_pt))
    for i_step in range(sample.n_load_step):
        txt_data             = np.loadtxt(ring.peak_dir+orient+'/'+sample.step_names[i_step]+'_peakfit_results.txt', dtype=float)
        l_cent[i_step, :], l_errs[i_step, :], l_amps[i_step, :]  = txt_data[:, 2], txt_data[:, 3], txt_data[:, 4]
        u_cent[i_step, :], u_errs[i_step, :], u_amps[i_step, :]  = txt_data[:, 5], txt_data[:, 6], txt_data[:, 7]
except:
    for i_step in range(sample.n_load_step):
        DA.write_fit_results(sample, ring, x1d, z1d, i_step, orient)
#%% 
        
##############################################################################
print('reading DIC..')
##############################################################################        

sigma_max = 0.05 

dic_list  = [ [[0], [0], [0]] ]
for i_step in range(1, sample.n_load_step):
    dic_path              = data_dir + sample.name + '/snapshots/dic/' + dic_files[i_step] + '.csv'
    dic_data              = DR.vic2d_reader(dic_path)
    good_correlation      = np.abs(dic_data[:, 4]) < sigma_max
    dic_x                 = dic_data[:, 5][good_correlation] - dic_center[0]
    dic_y                 = dic_data[:, 6][good_correlation] - dic_center[1]
    if orient == 'h':
        dic_strain            = dic_data[:,  9][good_correlation]
    if orient == 'v':
        dic_strain            = dic_data[:, 10][good_correlation]
    dic_list.append([dic_x, dic_y, dic_strain])

##############################################################################
print('filtering data..')
##############################################################################

# calculate ring radii
radii                 = (u_cent - l_cent) / 2  

# determine which points to use in analysis
err_max               = 0.5   # 2 norm of error / 2 norm of data
amp_min               = 200   # counts
use                   = np.ones((sample.n_load_step, sample.n_data_pt), dtype=bool)      
use[l_errs>err_max]   = False
use[u_errs>err_max]   = False
use[l_amps<amp_min]   = False
use[u_amps<amp_min]   = False  

# total variation filtering and dic data point matching
x, z, fits, t_strain  = [], [], [], []
xvals                 = np.unique(x1d)
for i_step in range(sample.n_load_step):
    s_x, s_z, s_fits, s_dic   = [], [], [], []
    for xval in xvals:
        use_col               =   use[i_step][x1d==xval] 
        x_col                 =           x1d[x1d==xval][use_col]
        z_col                 =           z1d[x1d==xval][use_col]
        radii_col             = radii[i_step][x1d==xval][use_col]
        if orient == 'h':
            fit_col               = DA.total_variation(radii_col, ring.lambda_h) 
        if orient == 'v':
            fit_col               = DA.total_variation(radii_col, ring.lambda_v) 
        dic_col               = DA.find_closest_vic2d(dic_list[i_step][0], dic_list[i_step][1], dic_list[i_step][2], x_col, z_col)
        s_x.append(x_col)
        s_z.append(z_col)
        s_fits.append(fit_col)
        s_dic.append(dic_col)
        path                  = ring.filt_dir+orient+'/'+'radii_fits_'+sample.step_names[i_step]+'_x_'+str(xval)+'.tiff'
        DA.plot_data_fit(path, radii_col, fit_col, sample.step_names[i_step]+', x='+str(xval))
    x.append(s_x)
    z.append(s_z)
    fits.append(s_fits)
    t_strain.append(s_dic)

##############################################################################
print('calculating strain..')
##############################################################################

#    find reference diameter    
refs                  = []
for i_step in range(sample.n_load_step):
    for i_col in range(len(fits[i_step])):
        refs.append(DA.find_sin_ref(fits[i_step][i_col], sample.detector_dist))
sin_ref              = np.mean(refs)

#    calculate elastic strains
e_strain              = []
for i_step in range(sample.n_load_step):
    s_e_strain = []
    for i_col in range(len(fits[i_step])):
        two_theta = np.arctan(fits[i_step][i_col] / sample.detector_dist)
        sin_theta = np.sin(two_theta / 2)
        s_e_strain.append( (sin_ref / sin_theta) - 1 )
    e_strain.append(s_e_strain)
    
#    calculate plastic strain
p_strain              = []
for i_step in range(sample.n_load_step):
    s_p_strain = []
    for i_col in range(len(fits[i_step])):
        s_p_strain.append(t_strain[i_step][i_col] - e_strain[i_step][i_col])
    p_strain.append(s_p_strain)
    
#    average strains
avg_e_strain = []
avg_p_strain = []
avg_t_strain = []
for i_step in range(sample.n_load_step):
    zvals                 = np.unique(np.concatenate(z[i_step]))
    row_avg_e_strain      = []
    row_avg_p_strain      = []
    row_avg_t_strain      = []
    for zval in zvals:
        row_e_strain          = []
        row_p_strain          = []
        row_t_strain          = []
        for i_col in range(len(e_strain[i_step])):
            if (z[i_step][i_col]==zval).any():
                row_e_strain.append(e_strain[i_step][i_col][z[i_step][i_col]==zval])
                row_p_strain.append(p_strain[i_step][i_col][z[i_step][i_col]==zval])
                row_t_strain.append(t_strain[i_step][i_col][z[i_step][i_col]==zval])
        row_avg_e_strain.append(np.mean(row_e_strain))
        row_avg_p_strain.append(np.mean(row_p_strain))
        row_avg_t_strain.append(np.mean(row_t_strain))
    avg_e_strain.append(row_avg_e_strain)
    avg_p_strain.append(row_avg_p_strain)
    avg_t_strain.append(row_avg_t_strain)

##############################################################################
print('plotting..')
##############################################################################

fig_size     = [6,4]         # width, height of figure in inches
legend_loc   = [1.4, 1.03]   # x, y coordinates of legend
label_size   = 18            # font size of x and y axis levels
tick_size    = 14            # font size of ticks
n_levels     = 21            # number of colorbar ticks to display
n_colorbar   = 11            # number of levels ot use in contour plot
scatter_size = 25            # size of scatter plot points
n_contour    = 1000          # number of points in each dimension to use in contour plot


def strain_plots(strain, avg_strain, strn_rng, descrip):

    # plot strain vs location for each column of data points
    for i_step in range(sample.n_load_step):
    
        plt.close('all')
        plt.figure(figsize=fig_size)
        plt.tick_params(labelsize=tick_size)
        plt.xlabel('y-coordinate', fontsize=label_size)
        plt.ylabel('strain',       fontsize=label_size)
        plt.ylim(strn_rng)
        for i_col in range(len(strain[i_step])):
            plt.plot(z[i_step][i_col], strain[i_step][i_col], '-', label='x = '+'%6.1f'%x[0][i_col][0], lw=2)
        plt.grid(True)
        plt.legend(numpoints=1, fontsize=tick_size, bbox_to_anchor=legend_loc)
        plt.savefig(ring.strn_dir+orient+'/'+descrip+'_strain_line_'+sample.step_names[i_step]+'.png', bbox_inches='tight', pad_inches=0.1)
    
    # plot strains averaged across columns
    plt.close('all')
    plt.figure(figsize=fig_size)
    plt.tick_params(labelsize=tick_size)
    plt.xlabel('y-coordinate', fontsize=label_size)
    plt.ylabel('strain',       fontsize=label_size)
    plt.ylim(strn_rng)
    for i_step in range(sample.n_load_step):
        plt.plot(np.unique(np.concatenate(z[i_step])), avg_strain[i_step], '-', label=sample.step_names[i_step]+' '*(10-len(sample.step_names[i_step])), lw=2)
    plt.grid(True)
    plt.legend(numpoints=1, fontsize=tick_size, bbox_to_anchor=legend_loc)
    plt.savefig(ring.strn_dir+orient+'/'+descrip+'_strain_line_averaged.png', bbox_inches='tight', pad_inches=0.1)

    # make 2d plots
    for i_step in range(sample.n_load_step):
        
        # contour and colorbar levels
        con_levels  = np.linspace(strn_rng[0], strn_rng[1], num=n_levels)
        cb_levels   = np.linspace(strn_rng[0], strn_rng[1], num=n_colorbar)
        
        # scatter plot
        plt.close('all')
        plt.figure(figsize=fig_size)
        plt.tick_params(labelsize=tick_size)
        plt.xlabel('x-coordinate', fontsize=label_size)
        plt.ylabel('y-coordinate', fontsize=label_size)
        for i_col in range(len(strain[i_step])):
            plt.scatter(x[i_step][i_col], z[i_step][i_col], s=scatter_size, c=strain[i_step][i_col], vmin=strn_rng[0], vmax=strn_rng[1], lw=0)
        cb = plt.colorbar(ticks=cb_levels, orientation='vertical')
        cb.ax.tick_params(labelsize=tick_size)
        plt.savefig(ring.strn_dir+orient+'/'+descrip+'_strain_scatter_'+sample.step_names[i_step]+'.png', bbox_inches='tight', pad_inches=0.1)
        
        # contour plot
        stepx, stepz, stepstrain     = [], [], []
        for i_col in range(len(strain[i_step])):
            for i_pt in range(len(strain[i_step][i_col])):
                stepx.append(x[i_step][i_col][i_pt])
                stepz.append(z[i_step][i_col][i_pt])
                stepstrain.append(strain[i_step][i_col][i_pt])
            
        xa               = np.linspace(np.min(stepx), np.max(stepx), num=n_contour, endpoint=True)
        za               = np.linspace(np.min(stepz), np.max(stepz), num=n_contour, endpoint=True)
        var_array        = griddata(stepx, stepz, stepstrain, xa, za, interp='linear')
        x_array, y_array = np.meshgrid(xa, za)    
        
        plt.close('all')
        plt.figure(figsize=fig_size)
        plt.tick_params(labelsize=tick_size)
        plt.xlabel('x-coordinate', fontsize=label_size)
        plt.ylabel('y-coordinate', fontsize=label_size)
        plt.contourf(x_array, y_array, var_array, con_levels, linestyle = 'solid')
        cb = plt.colorbar(ticks=cb_levels, orientation='vertical')
        cb.ax.tick_params(labelsize=tick_size)
        plt.savefig(ring.strn_dir+orient+'/'+descrip+'_strain_contour_'+sample.step_names[i_step]+'.png', bbox_inches='tight', pad_inches=0.1)
        plt.close('all')
        
strain_plots(e_strain, avg_e_strain, e_rng, 'elastic')
strain_plots(p_strain, avg_p_strain, p_rng, 'plastic')
strain_plots(t_strain, avg_t_strain, t_rng, 'total')
 
##############################################################################
print('writing..')
##############################################################################   

for i_step in range(sample.n_load_step):
    path = ring.strn_dir+orient+'/'+'strain_results_'+sample.step_names[i_step]+'.txt'
    out = open(path, 'w')
    for i_col in range(len(e_strain[i_step])):
        for i_data_pt in range(len(e_strain[i_step][i_col])):
            out.write('%22.16f'%x[i_step][i_col][i_data_pt]+'\t'+
                      '%22.16f'%z[i_step][i_col][i_data_pt]+'\t'+
                      '%22.16f'%e_strain[i_step][i_col][i_data_pt]+'\t'+
                      '%22.16f'%p_strain[i_step][i_col][i_data_pt]+'\t'+
                      '%22.16f'%t_strain[i_step][i_col][i_data_pt]+'\n')
    out.close()