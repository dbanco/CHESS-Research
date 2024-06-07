##################################################################################### SPOT FINDING
#####################################################################################
#####################################################################################

# Plots out the raw data and indicates where the spots are
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy.spatial.distance import cdist
from lmfit.models import Gaussian2dModel
from skimage.feature import blob_dog, blob_log, blob_doh
from scipy.spatial.distance import cdist
from skimage import io, segmentation, color
from skimage.measure import regionprops
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import lmfit
from lmfit.lineshapes import gaussian2d, lorentzian
from lmfit.models import Gaussian2dModel

def omegaToFrame(omega):
    step = 0.249788
    startPos = 0

    frame = round((math.degrees(omega) - startPos)/step) + 4
    return frame
    
def tthToRadius(tth):
    # stores detector to sample distance in pixels
    sampleDetector = 609.5449512495032
    mmPerPixel = 0.0748
    sampleDetector = sampleDetector/mmPerPixel
    radius = round(sampleDetector * math.tan(tth))
    return radius
    
def etaToPos(eta, radius, center):
        n, m = center
        n = n + round(radius * math.sin(eta))
        m = m + round(radius * math.cos(eta))
        return n, m

# gets all the angular coordinates of the spots in a grain and stores them into an array of [tth1, eta1, ome1, tth2, eta2, ...]
def getAngular(spots):
    posList = []
    for i in range(len(spots)):
        ID = spots.iloc[i, 0]
        ttheta = spots.iloc[i, 10]
        eta = spots.iloc[i, 11]
        ome = spots.iloc[i, 12]
        posList.append(int(ID))
        posList.append(float(ttheta))
        posList.append(float(eta))
        posList.append(float(ome))
    return posList

# returns a list with angular data values
def getData(rowCount, grains, grain):
    # range(rowCount) for all spots files
    angList = []
    
    zeroes = ""
    for j in range(5 - len(str(grain))):
        zeroes = zeroes + "0"
    spots = pd.read_csv('spots_' + zeroes + str(grain) + '.out', delim_whitespace = True)
    angList = getAngular(spots)
    return angList
        
# inputs (center of the image, eta angle at center of ring in degrees, ring radius in pixels)
# outputs image with marked spot
def locate(center, eta, tth):
    radius = tthToRadius(tth)
    position = etaToPos(eta, radius, center)
    n, m = position
    return position

# given a starting point, this code locates the nearest spot
def findBlob(x, y, image):
    input_point = (x, y)

    input_point = np.array([input_point])

    
    blobs_dog = blob_dog(image, min_sigma=1, max_sigma=30, threshold=200)
    
    if blobs_dog.size == 0:
        return 0, 0
    
    # obtain and store closest blob
    distances = cdist(input_point, blobs_dog[:, :2], 'euclidean')
    closest_blob_index = np.argmin(distances)
    closest_blob = blobs_dog[closest_blob_index]

    # Extract the (x, y) coordinates of the closest blob
    x, y, _ = closest_blob
    return x, y

# Fit a 2D Gaussian to the spot
def fitGaussian(image, x, y):
    # Create meshgrid for fitting
    y_vals, x_vals = np.indices(image.shape)

    # Define the Gaussian model
    gaussian_model = Gaussian2dModel()

    # Create parameters for the fit
    params = gaussian_model.guess(image.ravel(), x=x_vals.ravel(), y=y_vals.ravel())

    # Perform the fit
    result = gaussian_model.fit(image.ravel(), params, x=x_vals.ravel(), y=y_vals.ravel())
    return result

def gaussianFit(image):
    
    x, y = findBlob(image.shape[1]//2, image.shape[0]//2, image)
    print((x,y))
    result = fitGaussian(image, x, y)

    return result

def gaussianRecon(result,roi_size):
    fitted_params = result.params
    
    # Extract the parameters
    amplitude = fitted_params['amplitude'].value
    center_x = fitted_params['centerx'].value
    center_y = fitted_params['centery'].value
    sigma_x = fitted_params['sigmax'].value
    sigma_y = fitted_params['sigmay'].value
    
    recon = np.zeros((roi_size,roi_size))

    for x in range(roi_size):
        for y in range(roi_size):
            recon[y, x] = amplitude * np.exp(-((x - center_x)**2 / (2 * sigma_x**2) + (y - center_y)**2 / (2 * sigma_y**2)))
            
    return recon

def getGaussFitParams(result):
    fitted_params = result.params
    amplitude = fitted_params['amplitude'].value
    center_x = fitted_params['centerx'].value
    center_y = fitted_params['centery'].value
    sigma_x = fitted_params['sigmax'].value
    sigma_y = fitted_params['sigmay'].value
    
    factor = 2*np.sqrt(2*np.log(2))
    fwhm_x = factor*sigma_x
    fwhm_y = factor*sigma_y
    
    return np.array([amplitude,center_x,center_y,fwhm_x,fwhm_y])
    