"""
Data Reading Functions

@author: Kenny Swartz
06/07/2016
"""
import os, csv
import numpy as np


def vic2d_reader(path):
    """ function reads in a vic2d output .csv file
    
    inputs:
    path              : path to .csv file
    
    outputs:
    data              : numpy array of data (number of data points x columns of output data) """
    
    fid        = open(path, 'r')                         
    data       = []                                       
    for row in csv.reader(fid, delimiter=',', quotechar='|'): 
        if row != []:                               # skip blank rows
            data.append(row)                        
    fid.close()                                     
    
    data       = np.array(data[1:], dtype=float)          # skip header
    
    return data    



def ge2_reader(path, header_size=4096, image_size=2048):
    """ function reads in an image(s) from the ge2 detector 
    
    inputs:
    path              : path to image file 
    header_size       : size of image header in bits
    image_size        : size of each image in pixels 
    
    outputs:
    images            : numpy array of images (number of images x image size x image size) """
    
    fid        = open(path, 'r')                    
    image_1d   = np.fromfile(fid, dtype=np.uint16)  
    fid.close()                                     
    
    num_images = (image_1d.shape[0] - header_size) / image_size**2
    images     = np.array(image_1d[header_size:].reshape(num_images,image_size,image_size), dtype=float)
    
    return images
  
def ge2_reader_image(path,image_num, header_size=4096, image_size=2048):
    """ function reads in an image(s) from the ge2 detector 
    
    inputs:
    path              : path to image file 
    image_num         : index of image to read (starting at 0)
    header_size       : size of image header in bits
    image_size        : size of each image in pixels 
    
    outputs:
    images            : numpy array of images (number of images x image size x image size) """
    
    fid        = open(path, 'rb')
    fid.seek(header_size + 2*image_num*(image_size**2)  , os.SEEK_SET)                    
    image_1d   = np.fromfile(fid,count=image_size**2, dtype=np.uint16)  
    fid.close()                                     
    
    image     = np.array(image_1d.reshape(image_size,image_size), dtype=float)
    
    return image
  
def get_ge2_path(directory, dir_num, file_num, file_num_digits=5):
    """ function creates the path to a .ge2 file using Chess conventions
    
    inputs:
    directory         : directory where .ge2 image is located (1 above specific data point directory)
    dir_num           : specific data point directory number
    file_num          : number in .ge2 file name
    file_num_digits   : number of digits in .ge2 file name
    
    outputs:
    path              : path to desired .ge2 file """
    
    file_name  = 'ff_' + (file_num_digits-len(str(file_num)))*'0' + str(file_num) + '.ge2'
    path       = os.path.join(directory, str(dir_num),'ff', file_name)

    return path