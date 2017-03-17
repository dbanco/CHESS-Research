"""
Data Analysis Functions

@author: Kenny Swartz
06/07/2016
"""
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import os, importlib 
import DataReader, PeakFitting, Toolbox
reload(DataReader)
reload(PeakFitting)
reload(Toolbox)


class Specimen:              
    def __init__(self, name, data_dir, out_dir, step_names, dic_files, dark_dirs, init_dirs, dic_center, x_range, x_num, y_range, y_num, detector_dist, true_center, e_rng, p_rng, t_rng, E, G, v):      

        self.name             = name
        self.data_dir         = os.path.join(data_dir, name)
        self.out_dir          = os.path.join(out_dir, name)        
        
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)            
        
        self.step_names       = step_names
        self.n_load_step      = len(step_names)
        self.dic_files        = dic_files
        self.dark_dirs        = dark_dirs
        self.init_dirs        = init_dirs
        self.dic_center       = dic_center
        self.x_range          = x_range
        self.x_num            = x_num
        self.y_range          = y_range
        self.y_num            = y_num
        self.n_data_pt        = x_num*y_num
        self.detector_dist    = detector_dist
        self.true_center      = true_center
        self.e_rng            = e_rng
        self.p_rng            = p_rng
        self.t_rng            = t_rng
        self.E                = E
        self.G                = G
        self.v                = v

    def load_image(self,img_num,load_step):
        dir_start       = self.init_dirs[load_step]
        dir_num         = str(dir_start  + img_num)
        dir_num_dark    = str(self.dark_dirs[load_step])
        
        im_dir          = os.path.join(self.data_dir,dir_num,'ff')
        im_file         = os.listdir(im_dir)        
        assert len(im_file) == 1
        im_path = os.path.join(im_dir,im_file[0])
        
        dark_dir        = os.path.join(self.data_dir,dir_num_dark,'ff')
        dark_file       = os.listdir(dark_dir)
        assert len(dark_file) == 1
        dark_path       = os.path.join(dark_dir,dark_file[0])
    
        dark_image      = DataReader.ge2_reader_image(dark_path,0)
        if len(dark_image.shape) > 1:
            dark_image        = np.mean(dark_image, axis=0)        
        
        ring_image      = DataReader.ge2_reader_image(im_path,0) 
        
        img             = ring_image - dark_image
        return img        
        
class Ring:
    def __init__(self, sample, name, radius, dr, min_amp, vec_frac):
        
        self.name             = name
        self.ring_dir         = sample.out_dir+name+'\\'
        self.peak_dir         = self.ring_dir+'peak_fits\\'
        self.strn_dir         = self.ring_dir+'strain\\'
        
        if not os.path.exists(self.ring_dir):
            os.mkdir(self.ring_dir)
          
        if not os.path.exists(self.peak_dir):
            os.mkdir(self.peak_dir)       
            
        if not os.path.exists(self.strn_dir):
            os.mkdir(self.strn_dir) 
            
        self.radius           = radius
        self.dr               = dr
        self.min_amp          = min_amp
        self.vec_frac         = vec_frac
        

def strain_rosette_least_squares(two_theta_0, gamma, two_theta):
    """ function calculates exx, eyy, and exy from an overdetermined system of 2 theta measurements from diffraction image
    
    inputs:
    two_theta_0           : reference value correponding to the unstrained 2 theta value of a given ring
    gamma                 : array of angles corresponding to angle of the diffraction vector (after diffraction image is converted to r, gamma coordinates)
    two_theta             : array of two_theta peak fits corresponding to the gamma values above
    
    outputs:
    exx                   : normal strain component in x
    eyy                   : normal strain component in y
    exy                   : shear strain component in xy """
        
    # calculate normal strains
    normal_strains        = ( np.sin(two_theta_0 / 2) / np.sin( two_theta / 2) ) - 1
    
    # build least squares linear matrix to solve for 2D strain components
    A                     = np.zeros((gamma.shape[0], 3))
    A[:, 0]               = (1+np.cos(2*gamma)) / 2
    A[:, 1]               = (1-np.cos(2*gamma)) / 2
    A[:, 2]               = np.sin(2*gamma)
    
    # solve linear least squares problem
    exx, eyy, exy         = la.lstsq(A, normal_strains)[0]
    
    return exx, eyy, exy
    
    
def calculate_elastic_strain(sample, two_theta_0, gamma, two_theta, use, min_vecs):
    
    e_exx                 = np.zeros((sample.n_load_step, sample.n_data_pt))
    e_eyy                 = np.zeros((sample.n_load_step, sample.n_data_pt))
    e_exy                 = np.zeros((sample.n_load_step, sample.n_data_pt))
    
    for i_step in range(sample.n_load_step):
        for i_data_pt in range(sample.n_data_pt):
            use_gamma                = gamma    [i_step, i_data_pt][use[i_step, i_data_pt]]
            use_two_theta            = two_theta[i_step, i_data_pt][use[i_step, i_data_pt]]
            try:
                assert use_gamma.shape[0] >= min_vecs
                exx, eyy, exy            = strain_rosette_least_squares(two_theta_0, use_gamma, use_two_theta)
                e_exx[i_step, i_data_pt] = exx
                e_eyy[i_step, i_data_pt] = eyy
                e_exy[i_step, i_data_pt] = exy
            except:
                continue

    return e_exx, e_eyy, e_exy



def hookes_law(exx, eyy, exy, E, G, v):
    """ function calculates stress from strain using Hooke's law (assumes plane strain)
    
    inputs:
    exx                   : array of xx normal strain values 
    exy                   : array of xy shear strain  values 
    eyy                   : array of yy normal strain values 
    E                     : elastic modulus
    G                     : shear modulus
    v                     : poisson's ratio
    
    outputs:
    sxx                   : array of xx normal stress values
    syy                   : array of yy normal stress values
    sxy                   : array of xy shear stress values
    szz                   : array of zz normal stress values """
    
    # calculate stress
    sxx         = E * ( (1-v)*exx + v*eyy ) / ( (1+v)*(1-2*v) ) 
    syy         = E * ( (1-v)*eyy + v*exx ) / ( (1+v)*(1-2*v) ) 
    sxy         = G * exy
    szz         = v * (sxx + syy)
    
    return sxx, syy, sxy, szz



def write_2thetas(sample, ring, num_vecs, dgamma, x1d, y1d, step_num, fwhm0=10, amp0=500, plot_flag=False):

    # read in dark image
    dark_path             = sample.data_dir+str(sample.dark_dirs[step_num])+'\\ff\\'
    dark_file             = os.listdir(dark_path)
    assert len(dark_file) == 1
    dark_image            = DataReader.ge2_reader(dark_path+dark_file[0])
    if len(dark_image.shape) > 1:
        dark_image        = np.mean(dark_image, axis=0)

    # initialize storage arrays
    vec_gamma             = np.linspace(-np.pi+(dgamma/2), np.pi-(dgamma/2), num=num_vecs)
    two_theta             = np.zeros((sample.n_data_pt, num_vecs))
    peak_amps             = np.zeros((sample.n_data_pt, num_vecs))
    peak_errs             = np.zeros((sample.n_data_pt, num_vecs))
    
    # loop through each grid point on sample
    for i_data_pt in range(sample.n_data_pt):
       
       # read image
       dir_num               = sample.init_dirs[step_num] + i_data_pt
       path                  = sample.data_dir+str(dir_num)+'\\ff\\'
       file                  = os.listdir(path)
       assert len(file) == 1
       print('reading image ' + str(dir_num), 'x = '+str(x1d[i_data_pt]), 'y = '+str(y1d[i_data_pt])) 
       image                 = DataReader.ge2_reader(path+file[0])[0]  # only using first image because of shutter timing error
       image                -= dark_image                              # subtract dark image
      
       # generate coordinates of each pixel and calculate radius and vector angle
       x, y                  = np.meshgrid(np.arange(image.shape[1], dtype=float), np.arange(image.shape[0], dtype=float))
       x                    -= sample.true_center[1]
       y                    -= sample.true_center[0]
       radius                = np.sqrt( x**2 + y**2 )     # covert x,y coordinates into r,omega coordinates
       gamma                 = np.arctan2(y, x)           # covert x,y coordinates into r,omega coordinates
       
       # loop through each diffraction vector
       for i_vec in range(num_vecs):
            
            # grab slice of detector pixels that are within domega of desired omega
            img_slice                     =  image[np.abs(gamma-vec_gamma[i_vec]) < dgamma]
            r_slice                       = radius[np.abs(gamma-vec_gamma[i_vec]) < dgamma]
            
            # grab section of slice that is within dr of ring radius
            img_slice                     =  img_slice[np.abs(r_slice-ring.radius) < ring.dr]
            r_slice                       =    r_slice[np.abs(r_slice-ring.radius) < ring.dr]
            
            # sort selected pixels values by radial coordinate
            sorted_indices                = np.argsort(r_slice)
            sorted_r                      =   r_slice[sorted_indices]
            sorted_peak                   = img_slice[sorted_indices]
    
            # fit peak to sorted selected pixel values
            ctr_ind, lo_ind, hi_ind       = PeakFitting.get_peak_fit_indices(sorted_peak)
            peak_bg_rm, _                 = PeakFitting.RemoveBackground(sorted_r, sorted_peak, sorted_r[lo_ind], sorted_r[hi_ind])
            peak_fit, p_opt, err          = PeakFitting.fitPeak(sorted_r, peak_bg_rm, sorted_r[ctr_ind], fwhm0, amp0)
            
            # calculate 2 theta 
            opp                           = p_opt[0]
            adj                           = sample.detector_dist
            two_theta[i_data_pt, i_vec]   = np.arctan(opp/adj)
            
            # store peak amplitude and relative error
            peak_amps[i_data_pt, i_vec]   = p_opt[3]
            peak_errs[i_data_pt, i_vec]   = err
            
            if plot_flag:
                plt.close('all')
                fig = plt.figure()
                ax  = fig.add_subplot(111)
                ax.plot(sorted_r, sorted_peak, 'ok')
                ax.plot(sorted_r, peak_bg_rm,  'or')
                ax.plot(sorted_r, peak_fit,    '-r')
                ax.text(0.01, 0.92, 'ctr = '+str(opp), transform=ax.transAxes, color='k', fontsize=14)
                if err < 0.5:
                    ax.text(0.01, 0.85, 'err = '+str(err), transform=ax.transAxes, color='k', fontsize=14)
                else:
                    ax.text(0.01, 0.85, 'err = '+str(err), transform=ax.transAxes, color='r', fontsize=14)
                plt.savefig(ring.peak_dir+str(i_data_pt)+'_'+str(vec_gamma[i_vec])+'.png')
                plt.close('all')
                
       if plot_flag:
           plt.close('all')
           plt.imshow(image, vmin=0, vmax=200)
           plt.savefig(ring.peak_dir+str(i_data_pt)+'_image.png') 
           plt.close('all')

    # write data to a text file
    out_path = ring.peak_dir+sample.step_names[step_num]+'_peakfit_results.txt'
    out_file = open(out_path, 'w') 
    for i_data_pt in range(sample.n_data_pt):
        out_file.write('%24.16f'%x1d[i_data_pt]                + '\\t')
        out_file.write('%24.16f'%y1d[i_data_pt]                + '\\t')
        for i_vec in range(num_vecs):
            out_file.write('%24.16f'%vec_gamma[i_vec]             + '\\t')
            out_file.write('%24.16f'%two_theta[i_data_pt, i_vec]  + '\\t')
            out_file.write('%24.16f'%peak_amps[i_data_pt, i_vec]   + '\\t')
            out_file.write('%24.16f'%peak_errs[i_data_pt, i_vec]   + '\\t')
        out_file.write('\n')
    out_file.close()   
    

def find_closest_vic2d(dic_x, dic_y, dic_var, x, y):
    """ function finds the closest vic2d data point to desired data points described by x and y locations
    
    inputs:
    dic_x             : 1d array of x coordinates from dic analysis
    dic_y             : 1d array of y coordinates from dic analysis
    dic_var           : 1d array of desired variables from dic analysis
    x                 : 1d array of x coordinates of desired data points
    y                 : 1d array of y coordinates of desired data points

    outputs:
    var                : 1d array of desired variables at closest points to input x,y coordinates """
    
    if type(x) == np.ndarray:
        var               = np.zeros(len(x))
        
        for i_loc in range(x.shape[0]):
            dist_2            = (dic_x - x[i_loc])**2 + (dic_y - y[i_loc])**2
            var[i_loc]        = dic_var[np.argmin(dist_2)]    
            
    if type(x) == np.float64:
        dist_2            = (dic_x - x)**2 + (dic_y - y)**2
        var               = dic_var[np.argmin(dist_2)]
        
    return var
