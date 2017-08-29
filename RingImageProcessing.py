"""
Created on Sat Jun 11 16:08:40 2016

@author: Kenny Swartz, Dan Banco
"""
import numpy as np
import matplotlib.pyplot as plt

from numpy.linalg import norm 
from scipy.optimize import leastsq
from scipy.signal import argrelmax
from sklearn.linear_model import Lasso

import DataAnalysis as DA
import PeakFitting as peak
import EllipticModels as EM
import LassoSolvers as Lasso
import os
import time


def fit_circle_nonlin_lstsq(x, y, plot_flag=0):
    """
    Fits circle to data points
    inputs
            x           x-coordinates of points to fit
            y           y-coordinates of points to fit
            plot_flag   plots circle of true
            
    outputs
            xc          x-coordinate of center of circle
            yc          y-coordinate of center of circle
            r           radius of circle
    """

    def residual(params, x, y):
        xc, yc, r = params
        return (x-xc)**2 + (y-yc)**2 - r**2
    
    def Jacobian(params, x, y):           
        xc, yc, r = params
        J         = np.zeros((x.shape[0],params.shape[0]))
        J[:,0]    = -2*(x-xc)
        J[:,1]    = -2*(y-yc)
        J[:,2]    = -2*r
        return J
    
    # nonlinear least squares fit
    r_guess = np.sqrt(x[0]**2 + y[0]**2)/2
    params0 = np.array([x[0], y[0], r_guess])  # use first data point as initial guess
    params  = leastsq(residual, params0, args=(x,y), Dfun=Jacobian)[0]
    xc, yc, r = params    
    print('xc = '+'%11.8f'%xc, 'zc = '+'%11.8f'%yc, 'r = '+'%11.8f'%r)
   
    if(plot_flag):
        # plotting
        t         = np.linspace(0,2*np.pi,num=10000)
        xx        = xc + r*np.cos(t)
        yy        = yc + r*np.sin(t)
        plt.close('all')
        plt.plot(x,  y,  'ok', ms=10)         # input data points
        plt.plot(xx, yy, '-b', lw=2)          # nonlinear fit
        plt.plot(xc, yc, 'or', ms=10)         # nonlinear fit center
        plt.xlim([np.min(xx),np.max(xx)])
        plt.ylim([np.min(yy),np.max(yy)])
        plt.axes().set_aspect('equal')
        plt.xlabel('x')
        plt.ylabel('z')
        plt.show()
    
    return xc, yc, r

def fit_ellipse_lstsq(x,y):
    """
    Fits ellipse to data points
    inputs
            x           x-coordinates of points to fit
            y           y-coordinates of points to fit
            
    outputs
            a           x-axis radius parameter
            b           y-axis radius parameter
            orient      orientation angle of the ellipse
            X0          x-coordinate of center of ellipse before rotation
            Y0          y-coordinate of center of ellipse before rotation
            X0_in       x-coordinate of center of ellipse after rotation
            Y0_in       y-coordinate of center of ellipse after rotation

    """
    orientation_tol = 1e-3
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    x1 = x - mean_x
    y1 = y - mean_y
    
    X = np.array([x1**2, x1*y1, y1**2, x1, y1 ])
    C = np.sum(X,axis=1)
    A = np.dot(X,np.transpose(X))
    B = np.linalg.solve(A,C)
    
    a = B[0]
    b = B[1]
    c = B[2]
    d = B[3]
    e = B[4]
    
    if( min(abs(b/a),abs(b/c)) > orientation_tol ):
        orientation_rad = 0.5*np.arctan(b/(c-a))
        cos_phi = np.cos( orientation_rad )
        sin_phi = np.sin( orientation_rad )
        a = a*cos_phi**2 - b*cos_phi*sin_phi + c*sin_phi**2
        b = 0
        c = a*sin_phi**2 + b*cos_phi*sin_phi + c*cos_phi**2
        d = d*cos_phi - e*sin_phi
        e = d*sin_phi + e*cos_phi
        
        mean_x = cos_phi*mean_x - sin_phi*mean_y
        mean_y = sin_phi*mean_x + cos_phi*mean_y
    else:
        orientation_rad = 0
        cos_phi = 1
        sin_phi = 0
    test = a*c
    if( test < 0 ):
        status = 'hyperbola found'
        print( status + ': did not identify ellipse')
    elif( test == 0):
        status = 'parabola found'
        print( status + ': did not identify ellipse')
    elif( test > 0):
        status = ''
        if( a<0 ): 
            a = -a
            c = -c
            d = -d
            e = -e
            
        X0 = mean_x - d/2/a
        Y0 = mean_y - e/2/c
        F = 1 + (d**2)/(4*a) + (e**2)/(4*c)
        a = np.sqrt(F/a)
        b = np.sqrt(F/c)
#        long_axis = 2*max(a,b)
#        short_axis = 2*min(a,b)
        # Need to apply rotation to center point to get actual center
        X0_in = X0*cos_phi - Y0*sin_phi
        Y0_in = X0*sin_phi + Y0*cos_phi
        
        return [a,b,orientation_rad,X0,Y0,X0_in,Y0_in]
        
    else:
        print('Something went wrong')
        return 0

def fit_ellipse_at_center_nonlin_lstsq(x, y, xc, yc, plot_flag=0):
    """
    Fits ellipse to data points
    inputs
            x           x-coordinates of points to fit
            y           y-coordinates of points to fit
            xc          x-coordinate of center of ellipse
            yc          y-coordinate of center of ellipse
            plot_flag   plots circle of true
            
    outputs
            a           x-axis radius parameter
            b           y-axis radius parameter
    """
    def residual(params, x, y, xc, yc):
        a, b = params
        return b**2*(x-xc)**2 + a**2*(y-yc)**2 - (a*b)**2
    
    def Jacobian(params, x, y, xc, yc):           
        a, b = params
        J         = np.zeros((x.shape[0],params.shape[0]))
        J[:,0]    = -2*a*b**2 + a*(y-yc)**2
        J[:,1]    = -2*b*a**2 + b*(x-xc)**2
        return J
    
    # nonlinear least squares fit
    r_guess = np.sqrt(x[0]**2 + y[0]**2)/2
    params0 = np.array((r_guess,r_guess))  # use first data point as initial guess
    params  = leastsq(residual, params0, args=(x,y,xc,yc), Dfun=Jacobian)[0]
    a, b = params    
    print('a = '+'%11.8f'%a,   'b = '+'%11.8f'%b)
   
    if(plot_flag):
        # plotting
        t         = np.linspace(0,2*np.pi,num=10000)
        xx        = xc + a*np.cos(t)
        yy        = yc + b*np.sin(t)
        plt.close('all')
        plt.plot(x,  y,  'ok', ms=10)         # input data points
        plt.plot(xx, yy, '-b', lw=2)          # nonlinear fit
        plt.plot(xc, yc, 'or', ms=10)         # center
        plt.xlim([np.min(xx),np.max(xx)])
        plt.ylim([np.min(yy),np.max(yy)])
        plt.axes().set_aspect('equal')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
    
    return a, b
    

def fit_ellipse_nonlin_lstsq(x, y, plot_flag=0):
    """
    Fits ellipse to data points
    inputs
            x           x-coordinates of points to fit
            y           y-coordinates of points to fit
            plot_flag   plots circle of true
            
    outputs
            xc          x-coordinate of center of ellipse
            yc          y-coordinate of center of ellipse
            a           x-axis radius parameter
            b           y-axis radius parameter
    """
    def residual(params, x, y):
        xc, yc, a, b = params
        return b**2*(x-xc)**2 + a**2*(y-yc)**2 - (a*b)**2
    
    def Jacobian(params, x, y):           
        xc, yc, a, b = params
        J         = np.zeros((x.shape[0],params.shape[0]))
        J[:,0]    = -2*(x-xc)*b**2
        J[:,1]    = -2*(y-yc)*a**2
        J[:,2]    = -2*a*b**2 + a*(y-yc)**2
        J[:,3]    = -2*b*a**2 + b*(x-xc)**2
        return J
    
    # nonlinear least squares fit
    r_guess = np.sqrt(x[0]**2 + y[0]**2)/2
    params0 = np.array((x[0], y[0],r_guess,r_guess))  # use first data point as initial guess
    params  = leastsq(residual, params0, args=(x,y), Dfun=Jacobian)[0]
    xc, yc, a, b = params    
    print('xc = '+'%11.8f'%xc, 'zc = '+'%11.8f'%yc, 
          'a = '+'%11.8f'%a,   'b = '+'%11.8f'%b)
   
    if(plot_flag):
        # plotting
        t         = np.linspace(0,2*np.pi,num=10000)
        xx        = xc + a*np.cos(t)
        yy        = yc + b*np.sin(t)
        plt.close('all')
        plt.plot(x,  y,  'ok', ms=10)         # input data points
        plt.plot(xx, yy, '-b', lw=2)          # nonlinear fit
        plt.plot(xc, yc, 'or', ms=10)         # nonlinear fit center
        plt.xlim([np.min(xx),np.max(xx)])
        plt.ylim([np.min(yy),np.max(yy)])
        plt.axes().set_aspect('equal')
        plt.xlabel('x')
        plt.ylabel('z')
        plt.show()
    
    return a, b, xc, yc

def n_std_threshold(img,n):
    """
    Thresholds image with threshold set to mean + n*standard_deviation
    inputs: img           image
            n             constant (threshold = mu + n*sigma)
            
    outputs: thresh_img   binary image    
    """
    mu = np.mean(img.flatten())
    sigma = np.std(img.flatten())
    
    thresh_low = mu + n*sigma

    return img > thresh_low
    

def distance_threshold(img,radius_low,radius_high,center=0):
    """ 
    Creates binary image by thresholding pixel distance from given center
    inputs
            img             image to be thresholded by pixel distance
            radius_low      distance above which pixels are kept 
            radius_high     distance below which pixels are kept
            center          [row,col] position from which distance is measured
    outputs
            bin_img         binary image where pixels equal to 1 satisfy
                            radius_low < distance < radius_high
    """
    bin_img = np.zeros(img.shape)
    n,m = img.shape
    if(center == 0):
        center = [n/2.0,m/2.0] 
        
    # limit to single ring
    for i in range(n):
        for j in range(m):
            pos = (i-center[0])**2 + (j-center[1])**2
            if ((pos > (radius_low)**2) & (pos < (radius_high)**2)):
                bin_img[i,j] = 1                         
    
    return bin_img == 1


def get_points(img,bin_img):
    """
    Gets coordinates and values of img from true values in binary image
    inputs
            img         an image
            bin_img     a binary image
            
    outputs
            x           x coordinate of pixels
            y           y coordinate of pixels
            f           value of pixels
    """
    y,x = np.nonzero(bin_img)
    f = img[y,x]
    return x,y,f


def radial_projection(img,center,r,num_r,theta,r_in,r_out,pf1=0):
    """
    Interpolates along a line in the radial direction
    
    inputs
            img     image
            center  center of the rings in the image
            r       radius of ring of interest
            num_r   number of points to sample along line
            theta   angle of line
            r_in    inside radius containing ring
            r_out   outside radius containing ring
            
    outputs
            r_project  image values at points along line
            r_domain   domain over which image values are defined
    """
    n,m = img.shape
    if(center==0): center = [round(n/2.0),round(m/2.0)]  

    r_domain  =   np.linspace(r_in,r_out,num_r)
    x_domain = r_domain*np.cos(theta) + center[0]
    y_domain = r_domain*np.sin(theta) + center[1]
    
    if(pf1):
        plt.figure(pf1)
        plt.plot(x_domain,y_domain,'-ow')
    
    r_project = np.zeros(r_domain.shape)      
    
    for ridx in range(len(r_domain)):
        
        # identify surrounding four points
        x = x_domain[ridx]
        y = y_domain[ridx]     
        x1 = np.floor( x )
        x2 = np.ceil(  x )
        y1 = np.floor( y )
        y2 = np.ceil(  y )
  
        # make sure we are in image
        if( (x1 < n) & (x2 < n) & (x1 > 0) & (x2 > 0) &
            (y1 < m) & (y2 < m) & (y1 > 0) & (y2 > 0) ):
            if( ((x2-x1) == 0) & ((y2-y1) == 0) ):
                r_project[ridx] = img[x1,y1]
            elif(x2-x1 == 0):
                r_project[ridx] = img[x1,y1] + \
                                (img[x2,y2]-img[x1,y1])*(y-y1)/(y2-y1)
            elif(y2-y1 == 0):
                r_project[ridx] = img[x1,y1] + \
                                (img[x2,y2]-img[x1,y1])*(x-x1)/(x2-x1)
            else:
                
                # interpolate
                a = np.matrix([x2-x,x-x1])
                Q = np.matrix([[img[x1,y1],img[x1,y2]],
                              [img[x2,y1],img[x2,y2]]])      
                b = np.matrix([[y2-y],[y-y1]])
                r_project[ridx] = np.dot(np.dot(a,Q),b)/((x2-x1)*(y2-y1))

    return r_project, r_domain
    
def line_normal_to_ellipse(x0,y0,a,b,xc,yc,distance):
    """
    Finds two points along line perpendicular to ellipse a distance from the
    point on the ellipse at angle theta about the center. Returns slope and 
    intercept of the line as well
    inputs
        x0          x-coordinate of points
        y0          y-coordinate of points
        distance    distance(s) from point to find two points along line
    outputs
        x           coordinates of two points
        y           coordinates of two points
    """
    x0 -= xc
    y0 -= yc
    if( np.abs(y0) > np.abs(x0)):
        # Slope of line normal to ellipse at that point dy/dx
        slope  = (y0*a**2)/(x0*b**2)
 
        if(np.isinf(slope**2)):
            y = np.array([y0-distance,y0+distance])
            x = np.array([x0,x0])
        else:              
            val = np.sqrt(distance**2/(1+slope**2))
            x =  np.array([x0 + val, x0 - val])
            y = slope*(x - x0) + y0
    else:
        # Slope of line normal to ellipse at that point dx/dy
        slope  = (x0*b**2)/(y0*a**2)

        if(np.isinf(slope**2)):
            x = np.array([x0-distance,x0+distance])
            y = np.array([y0,y0])
        else:
            val = np.sqrt(distance**2/(1+slope**2))
            y = np.array([y0 - val, y0 + val])
            x = slope*(y - y0) + x0
    return x + xc, y + yc
    

def radial_projection_ellipse(img,xc,yc,a,b,theta,distance,pf1):
    """
    Interpolates along a line in the radial direction
    
    inputs
            img      image
            center   center of the rings in the image
            a,b      major/minor axis of ring of interest
            theta    angle of line
            distance range over which to project radially
            d_in     inside window containing ring
            d_out    outside window containing ring]
            pf1      plots all interpolated points if true
    outputs
            f     image values at points along line
            x,y   domain over which image values are defined
    """
    # x-y coordinates of point on ellipse
    x0 = a*np.cos(theta) + xc
    y0 = b*np.sin(theta) + yc

    # Get endpoints of line perpendicular to ellipse at x0, y0
    x_l, y_l = line_normal_to_ellipse(x0,y0,a,b,xc,yc,distance)
    
    if(y_l[0] > 5000):print(theta)
    sort = (x_l-xc)**2 + (y_l-yc)**2
    sort_ind = [np.argmin(sort), np.argmax(sort)]
    
    x_line = x_l[sort_ind]
    y_line = y_l[sort_ind]                       

    
    if( x_line[1] == x_line[0]):
        x_domain = x_line[0]*np.ones(2*distance)
        y_domain = np.linspace(y_line[0],y_line[1],2*distance)
    elif( y_line[1] == y_line[0]):
        y_domain = y_line[0]*np.ones(2*distance)
        x_domain = np.linspace(x_line[0],x_line[1],2*distance)
    else:
        x_domain = np.linspace(x_line[0],x_line[1],2*distance)
        y_domain = np.linspace(y_line[0],y_line[1],2*distance)
        
    f = np.zeros(x_domain.shape)

    if(pf1):
        plt.figure(pf1)
        plt.plot(x_domain,y_domain,'-ow')

    for idx in range(len(x_domain)):
        
        # identify surrounding four points
        x = x_domain[idx]
        y = y_domain[idx]    
        x1 = np.floor( x )
        x2 = np.ceil(  x )
        y1 = np.floor( y )
        y2 = np.ceil(  y )
  
        # make sure we are in image
        n,m = img.shape
        if( (x1 < n) & (x2 < n) & (x1 > 0) & (x2 > 0) &
            (y1 < m) & (y2 < m) & (y1 > 0) & (y2 > 0) ):
            if( ((x2-x1) == 0) & ((y2-y1) == 0) ):
                f[idx] = img[y1,x1]
            elif(x2-x1 == 0):
                f[idx] = img[y1,x1] + (img[y2,x2]-img[y1,x1])*(y-y1)/(y2-y1)
            elif(y2-y1 == 0):
                f[idx] = img[y1,x1] + (img[y2,x2]-img[y1,x1])*(x-x1)/(x2-x1)
            else: 
                # interpolate
                a = np.matrix([x2-x,x-x1])
                Q = np.matrix([[img[y1,x1],img[y1,x2]],
                              [img[y2,x1],img[y2,x2]]])      
                b = np.matrix([[y2-y],[y-y1]])
                f[idx] = np.dot(np.dot(a,Q),b)/((x2-x1)*(y2-y1))

    return f, x_domain, y_domain

def azimuthal_projection(img,center,r,theta_1,theta_2,num_theta):
    """
    Interpolates along a line in the radial direction
    
    inputs
            img         image
            center      center of the rings in the image
            r           radius of ring of interest     
            num_theta   number of samples along theta
            
    outputs
            theta_project  image values at points along line
            theta_domain   domain over which image values are defined
    """
    n,m = img.shape
    if(center==0): center = [round(n/2.0),round(m/2.0)]  

    theta_domain  =   np.linspace(theta_1,theta_2,num_theta)
    theta_project = np.zeros(theta_domain.shape)
    
    for tidx in range(len(theta_domain)):
        # identify surrounding four points
        x = r*np.cos(theta_domain[tidx]) + center[0]
        y = r*np.sin(theta_domain[tidx]) + center[1]     
        x1 = int(np.floor( x ))
        x2 = int(np.ceil(  x ))
        y1 = int(np.floor( y ))
        y2 = int(np.ceil(  y ))
  
        # make sure we are in image
        if( (x1 < n) & (x2 < n) & (x1 > 0) & (x2 > 0) &
            (y1 < m) & (y2 < m) & (y1 > 0) & (y2 > 0) ):
            if((x2-x1 == 0) & (y2-y1 == 0)):
                theta_project[tidx] = img[y1,x1]
            elif((x2-x1) == 0):
                theta_project[tidx] = img[y1,x1] + \
                                (img[y2,x2]-img[y1,x1])*(y-y1)/(y2-y1)
            elif((y2-y1) == 0):
                theta_project[tidx] = img[x1,y1] + \
                                (img[y2,x2]-img[y1,x1])*(x-x1)/(x2-x1)
            else:
                
                # interpolate
                a = np.matrix([x2-x,x-x1])
                Q = np.matrix([[img[y1,x1],img[y1,x2]],
                              [img[y2,x1],img[y2,x2]]])      
                b = np.matrix([[y2-y],[y-y1]])
                theta_project[tidx] = np.dot(np.dot(a,Q),b)/((x2-x1)*(y2-y1))

    return theta_project, theta_domain


def azimuthal_projection_ellipse(img,center,a,b,rot,theta_1,theta_2,num_theta):
    """
    Interpolates along an ellipse in the azimuthal direction
    
    inputs
            img         image
            center      center of the rings in the image
            a           x-axis radius parameter of ellipse
            b           y-axis radius parameter of ellipse     
            rot         orientation of ellipse
            theta_1     first angle at which ellipse is sampled
            theta_2     last angle at which ellipse is sample
                        Angular distance between theta1 and theta2 should be
                        dtheta ot acheive uniformly spaced points around the 
                        entire ellipse
                        ie: 0 and 2*np.pi-dtheta
            num_theta   number of samples along theta
            
    outputs
            theta_project  image values at points along line
            theta_domain   domain over which image values are defined
    """
    n,m = img.shape
    if(center==0): center = [round(n/2.0),round(m/2.0)]  

    theta_domain  =   np.linspace(theta_1,theta_2,num_theta)
    theta_project = np.zeros(theta_domain.shape)
    
    xa = a*np.cos(theta_domain) + center[0]
    ya = b*np.sin(theta_domain) + center[1]    
    x  =  (xa)*np.cos(rot) + (ya)*np.sin(rot)
    y  = -(xa)*np.sin(rot) + (ya)*np.cos(rot) 
    
    for tidx in range(len(theta_domain)):
        # identify surrounding four points
        x_tidx = x[tidx]
        y_tidx = y[tidx]
        x1 = np.floor( x_tidx )
        x2 = np.ceil(  x_tidx )
        y1 = np.floor( y_tidx )
        y2 = np.ceil(  y_tidx )

        # make sure we are in image
        if( (x1 < n) & (x2 < n) & (x1 > 0) & (x2 > 0) &
            (y1 < m) & (y2 < m) & (y1 > 0) & (y2 > 0) ):
            if((x2-x1 == 0) & (y2-y1 == 0)):
                theta_project[tidx] = img[x1,y1]
            elif((x2-x1) == 0):
                theta_project[tidx] = img[x1,y1] + \
                                (img[y2,x2]-img[y1,x1])*(y_tidx-y1)/(y2-y1)
            elif((y2-y1) == 0):
                theta_project[tidx] = img[x1,y1] + \
                                (img[y2,x2]-img[y1,x1])*(x_tidx-x1)/(x2-x1)
            else:
                
                # interpolate
                a = np.matrix([x2-x_tidx,x_tidx-x1])
                Q = np.matrix([[img[y1,x1],img[y1,x2]],
                              [img[y2,x1],img[y2,x2]]])      
                b = np.matrix([[y2-y_tidx],[y_tidx-y1]])
                theta_project[tidx] = np.dot(np.dot(a,Q),b)/((x2-x1)*(y2-y1))

    return theta_project, theta_domain  
    
def azimuthal_projection_ellipse_radial_sum(img,center,a,b,theta_domain,radial_distance,pf1=0):
    """
    Interpolates along radial lines normal to an ellipse that are evenly spaced
    in azimuthal angle around the ellipse. Radial lines are averaged to
    compute an intensity value at each azimuthal angle along the perimeter of
    the ellipse. The exact sample points can be provided.
    
    inputs
            img            image
            center         center of the rings in the image
            a                  x-axis radius parameter of ellipse
            b                  y-axis radius parameter of ellipse     
            theta_domain       azimuthal domain where image values are 
                               interpolated
            radial_distance    distance from perimeter of ellipse over which
                               pixel values are averaged
            
    outputs
            theta_project  image values at points along line
            
    """
    n,m = img.shape
    if(center==0): center = [round(n/2.0),round(m/2.0)]  
    
    num_theta = len(theta_domain)
    f_az = np.zeros(theta_domain.shape)
    
    if(pf1):
        plt.figure(pf1)
        xx, zz = gen_ellipse(a,b,center[0],center[1],0,theta_domain)
        plt.plot(xx,zz,'w-*',markersize=10)
        plt.xlim([np.min(xx)-20, np.max(xx)+20])
        plt.ylim([np.min(zz)-20, np.max(zz)+20])
    
    for az_i, az in enumerate(theta_domain):
        f_rad, x_rad, y_rad = radial_projection_ellipse(img,
                                                        center[0],center[1],
                                                        a,b,az,
                                                        radial_distance,pf1)                                                    
        f_az[az_i] = np.mean(f_rad)

    return f_az, theta_domain
    
def gaussian_convolution(signal,sigma,C=4):
    """
    inputs
            signal  a one dimensional signal (n,) numpy array
            sigma   a standard deviation
            C       constant to determine length of gaussian kernel
            
    outputs
            filtered signal     convolution of gaussian kernel and input signal
                                where signal length is preserved
    """
    
    M = int(round(C*sigma + 1))
    
    gaussian = np.zeros(2*M+1)
    amp = 1/(np.sqrt(2*np.pi)*sigma)
    const = 1/(2*sigma)
    
    for m in range(-M,M+1):
        gaussian[m+M] = np.exp(-const*m**2)
    
    G = amp*gaussian
    
    return np.convolve(signal,G,mode='same')
    

def find_scale_space_maxima(signal,sigma,C=4,octaves=4):
    """
    inputs
            signal  a one dimensional signal (n,) numpy array
            sigma   a standard deviation
            C       constant to determine length of gaussian kernel
            octaves number of scales over which to find maxima
            
    outputs
            maxima  list of arrays containing locations of maxima at each scale
    """
    maxima = []
    maxima.append(argrelmax(signal))
    plt.figure(0)
    plt.plot(signal)
    for i in range(octaves):
       filt_signal = gaussian_convolution(signal,2*(i+1)*sigma,C)
       maxima.append(argrelmax(filt_signal)[0])
#       plt.figure(i+1)
       plt.plot(filt_signal)
    return maxima
        
def do_peak_fit(y,edgeCut,index=0,param=0,plot_flag=0):
    """
    Fits peak with two-sided Gaussian. Assumes peak is in the center and that
    edgeCut tells which portion of y to treat as background. interpolates
    between the points at the cutoff to subtract background before fitting
    a peak.
    
    inputs
            y           a 1D vector containing a peak
            edgeCut     fraction at edge of y to treat as background
            index       used to identify plot, appears in title
            param       used to identify plot, appears in title
            plot_flag   produces and saves plot to "plots" directory
    outputs
            mu          mean of the fitted two-sided Gaussian
            variance    total variance of the two-sided Gaussian
            amp         amplitude of the two-sided Gaussian
            err         ||fit - y||/||y||
    """
    amp_est                 = np.max(y)
    fwhm_est                = len(y)/2.0
    fit_domain              = np.arange(len(y))
    peakCtr, loCut, hiCut   = DA.get_peak_fit_indices(y,lo=edgeCut,hi=1-edgeCut)
    data_rm_back, back      = peak.RemoveBackground(fit_domain,y,loCut,hiCut)
    _, param_opt, err       = peak.fitPeak(fit_domain, data_rm_back, peakCtr, 
                                       fwhm=fwhm_est, amp=amp_est,
                                       FitType='Gaussian', n=1)
                                       
    mu, sigmaL, sigmaR, amp = param_opt
    
    if(plot_flag):
        plt.close(1)
        plt.figure(1)
        plt.plot(y,'o-b') 
        x_fit = np.linspace(fit_domain[0],fit_domain[-1],500)
        plt.plot(x_fit,peak.gaussian(x_fit,param_opt)+np.mean(back),'-r')   
        title = str(index) + '_err_' + str(err) + '_param_' + str(param)
        plt.title(title)
        plt.savefig(os.path.join('plots',title + '.png'))

    
    return mu, ((sigmaL+sigmaR)/2)**2, amp, err
    
def crop_around_max(x):
    """
    Finds max of vector and crops segment so max is at center
    inputs
            x           a 1D numpy array
    outputs
            x_cropped   cropped 1D numpy array
            indices     indices such that x_cropped = x[indices]
            center      index where center is located in x
    """
    out = np.argmax(x)
    center = round(out)
    if( center < (len(x)-1-center) ):
        x_crop = x[0:(2*center+2)]
        indices = range(0,int(2*center+2))
    elif( center > (len(x)-1-center) ):
        x_crop = x[(2*center-(len(x)-1)):len(x)]
        indices = range(int(2*center-(len(x)-1)),len(x))
    else:
        x_crop = np.copy(x)
        indices = range(0,len(x))
     
    return x_crop, indices, center
        
def ring_fit_nonlin(img,r_est,dr,center=0,n_std=3,pf1=False):
    """  
    Fits an ellipse to a ring of an xray diffraction image 
    inputs
        img-    image of xray diffraction rings    
        r_est-  estimate of radius of ring of interest in pixels
        dr-     space around ring on either side in pixels
        center          [row,col] position from which distance is measured
    outputs
        a-      x-axis radius parameter of fitted ellipse     
        b-      z-axis radius parameter of fitted ellipse
        xc-     x-coordinate of center of ellipse
        zc-     z-coordinate of center of ellipse
        rot-    orientation angle of ellipse
        
        
    EX:
    # Estimate radius of ring
    plt.close('all')
    plt.figure(1)
    plt.imshow(img,cmap='jet',vmin=0,vmax=200)
    xc = 1026
    zc = 1018
    radius_test = 606
    plt.plot([xc,xc+radius_test],[zc,zc],'o-w')
    
    # Fit ring
    radius = 370
    dr     = 30
    a,b,xc,zc,rot = RingIP.ring_fit(img,radius,dr)

    """  
    if(center == 0):
        n,m = img.shape
        center = [n/2.0,m/2.0] 
    if(pf1):
        plt.close(pf1)
        plt.figure(pf1)
        plt.imshow(img,cmap='jet',vmin=0,vmax=200,interpolation='nearest')
    
    r_in = r_est - dr
    r_out = r_est + dr

    # Keeps only pixels near ring
    bin_img = distance_threshold(img,r_in,r_out,center)    
    
    # Thresholds remainder of image with threshold 
    # 2 standard deviations above the mean
    thresh_img = n_std_threshold(img*bin_img,n_std)
    
    if(pf1):
        plt.close(pf1+1)
        plt.figure(pf1+1)
        plt.imshow(img*bin_img,cmap='jet',vmin=0,vmax=200,interpolation='nearest')    
    
        plt.close(pf1+2)
        plt.figure(pf1+2)
        plt.imshow(img*bin_img*thresh_img,cmap='jet',vmin=0,vmax=200,interpolation='nearest')  
        
        plt.close(pf1+3)
        plt.figure(pf1+3)
        plt.imshow(img*bin_img,cmap='jet',vmin=0,vmax=200,interpolation='nearest')    
    # Converts identified pixels to a set of x,z,f data points
    x,z,f = get_points(img,bin_img*thresh_img)

    # Fit ellipse to found points (ignores intensity)
    out = fit_ellipse_nonlin_lstsq(x,z)
    
    a = out[0]
    b = out[1]
    xc  = out[2]
    zc = out[3]

    if(pf1):
        plt.figure(pf1+3)
        xx,zz = gen_ellipse(a,b,xc,zc,0)
        plt.plot(xx,zz,'w-')
        plt.xlim([np.min(x)-10, np.max(x)+10])
        plt.ylim([np.min(z)-10, np.max(z)+10])

    return a, b, xc, zc     
    
    
def ring_fit(img,r_est,dr,center=0,n_std=2,pf1=False):
    """  
    Fits an ellipse to a ring of an xray diffraction image 
    inputs
        img-    image of xray diffraction rings    
        r_est-  estimate of radius of ring of interest in pixels
        dr-     space around ring on either side in pixels
        center          [row,col] position from which distance is measured
    outputs
        a-      x-axis radius parameter of fitted ellipse     
        b-      z-axis radius parameter of fitted ellipse
        xc-     x-coordinate of center of ellipse
        zc-     z-coordinate of center of ellipse
        rot-    orientation angle of ellipse
        
        
    EX:
    # Estimate radius of ring
    plt.close('all')
    plt.figure(1)
    plt.imshow(img,cmap='jet',vmin=0,vmax=200)
    xc = 1026
    zc = 1018
    radius_test = 606
    plt.plot([xc,xc+radius_test],[zc,zc],'o-w')
    
    # Fit ring
    radius = 370
    dr     = 30
    a,b,xc,zc,rot = RingIP.ring_fit(img,radius,dr)

    """  
    if(center == 0):
        n,m = img.shape
        center = [n/2.0,m/2.0] 
    if(pf1):
        plt.close(pf1)
        plt.figure(pf1)
        plt.imshow(img,cmap='jet',vmin=0,vmax=200)
    
    r_in = r_est - dr
    r_out = r_est + dr

    # Keeps only pixels near ring
    bin_img = distance_threshold(img,r_in,r_out,center)    
    
    # Thresholds remainder of image with threshold 
    # 2 standard deviations above the mean
    thresh_img = n_std_threshold(img*bin_img,n_std)

    # Converts identified pixels to a set of x,z,f data points
    x,z,f = get_points(img,bin_img*thresh_img)
    print('Number of points to fit: ' + str(len(x)))
    # Fit ellipse to found points (ignores intensity)
    out = fit_ellipse_lstsq(x,z)
    a = out[0]
    b = out[1]
    rot = out[2]
    xc = out[3]
    zc = out[4]

    if(pf1):
        plt.figure(pf1)
        plt.plot(x,z,'wx')
        xx,zz = gen_ellipse(a,b,xc,zc,rot)
        plt.plot(xx,zz,'w-')
        plt.xlim([np.min(x)-10, np.max(x)+10])
        plt.ylim([np.min(z)-10, np.max(z)+10])

    return a, b, xc, zc, rot     
     
def two_stage_ring_fit(img,edgeCut,r_est,dr,pf1=False,pf2=False):
    """  
    Fits an ellipse to a ring of an xray diffraction image 
    inputs
        img-    image of xray diffraction rings    
        edgeCut-edge cutoff for peak fitting (fraction of vector length)
        r_est-  estimate of radius of ring of interest in pixels
        dr-     space around ring on either side in pixels
        pf1-    plots output of 2nd stage in figure(pf1)
        pf2-    plots output of 1st stage in figure(pf2)
    outputs
        a-      x-axis radius parameter of fitted ellipse     
        b-      z-axis radius parameter of fitted ellipse
        xc-     x-coordinate of center of ellipse
        zc-     z-coordinate of center of ellipse
        rot-    orientation angle of ellipse
        
        
    EX:
    # Estimate radius of ring
    plt.close('all')
    plt.figure(1)
    plt.imshow(img,cmap='jet',vmin=0,vmax=200)
    xc = 1026
    zc = 1018
    radius_test = 606
    plt.plot([xc,xc+radius_test],[zc,zc],'o-w')
    
    # Fit ring
    edgeCut = 0.2
    radius = 370
    dr     = 30
    a,b,xc,zc,rot = RingIP.two_stage_ring_fit(img,edgeCut,radius,dr)

    """  
    if(pf1):
        plt.close(pf1)
        plt.figure(pf1)
        plt.imshow(img,cmap='jet',vmin=0,vmax=200)
        
    if(pf2):
        plt.close(pf2)
        plt.figure(pf2)
        plt.imshow(img,cmap='jet',vmin=0,vmax=200)
    
    r_in = r_est - dr
    r_out = r_est + dr

    # Keeps only pixels near ring
    bin_img = distance_threshold(img,r_in,r_out)    
    
    # Thresholds remainder of image with threshold 
    # 2 standard deviations above the mean
    thresh_img = n_std_threshold(img*bin_img,1)

    # Converts identified pixels to a set of x,z,f data points
    x,z,f = get_points(img,bin_img*thresh_img)

    # Fit ellipse to found points (ignores intensity)
    out1 = fit_ellipse_lstsq(x,z)
    a1 = out1[0]
    b1 = out1[1]
    rot1 = out1[2]
    xc1 = out1[3]
    zc1 = out1[4]

    if(pf2):
        plt.figure(pf2)
        plt.plot(x,z,'wx')
        xx,zz = gen_ellipse(a1,b1,xc1,zc1,rot1)
        plt.plot(xx,zz,'w-')
        plt.xlim([np.min(xx)-10, np.max(xx)+10])
        plt.ylim([np.min(zz)-10, np.max(zz)+10])

    # Identify approximate number of significant samples to interpolate
    num_r = round(r_out-r_in)
    # Identify smallest possible angular interval dtheta
    r1 = np.mean(np.array([a1,b1]))
    dtheta = np.abs(np.pi/4 - np.arctan((r1+1)/(r1-1)))
    # Might include integer number of steps as an input parameter
    step = 20*dtheta
    theta = np.arange(dtheta,2*np.pi + step,step)
    
    mu = np.zeros(theta.shape)
    err = np.zeros(theta.shape)
    x2 = np.array([])
    z2 = np.array([])
    for i in range(len(theta)-1):
        # Get points along line normal to ellipse at angle theta[i]
        f, x_dom, y_dom = radial_projection_ellipse(img,xc1,zc1, 
                                                    a1,b1,rot1,num_r,
                                                    [theta[i-1],theta[i],
                                                    theta[i+1]],dr)
        # Crop segment so peak is at the center                                                    
        f_crop, crop_ind, _ = crop_around_max(f)
        
        x_dom_crop = x_dom[crop_ind]
        y_dom_crop = y_dom[crop_ind]
        
        # Fit peaks in cropped segments if cropped segment is long enough
        if (len(f_crop) > 0.2*len(f)):
            mu[i],_,_,err[i] = do_peak_fit(f_crop,edgeCut)
        else:
            mu[i] = 0
            err[i] = 1
        
        # Keep only peak fits within error threshold
        if(err[i]<0.2):
            xmid = x_dom_crop[int(round(mu[i]))]
            ymid = y_dom_crop[int(round(mu[i]))]
            x2 = np.append(x2,xmid)
            z2 = np.append(z2,ymid)
            
            if(pf1):
                x1p  = x_dom_crop[0]
                x2p  = x_dom_crop[-1]
                y1p  = y_dom_crop[0]
                y2p  = y_dom_crop[-1]
                plt.figure(pf1)
                plt.plot((x1p,xmid,x2p),(y1p,ymid,y2p),'o-w')

    
    # Fit ellipse to the top of the peaks
    out2 = fit_ellipse_lstsq(x2,z2)
    a2 = out2[0]
    b2 = out2[1]
    rot2 = out2[2]
    xc2 = out2[3]
    zc2 = out2[4]

    if(pf1):
        plt.figure(pf1)
        xx,zz = gen_ellipse(a1,b1,xc1,zc1,rot1)
        plt.plot(xx,zz,'w-')
        plt.xlim([np.min(xx)-10, np.max(xx)+10])
        plt.ylim([np.min(zz)-10, np.max(zz)+10])
    return a2, b2, xc2, zc2, rot2

def gen_ellipse(a,b,xc,yc,rot,domain=0):
    """
    Generates x,y data points along the perimeter of ellipse defined by inputs
    inputs            
            a-      x-axis radius parameter of fitted ellipse     
            b-      z-axis radius parameter of fitted ellipse
            xc-     x-coordinate of center of ellipse
            yc-     y-coordinate of center of ellipse
            domain- theta domain at which to plot points on ellipse
    
    outputs
            x       x-coordinates of points on ellipse
            y       y-coordinates of points on ellipse
    """
    if(np.isscalar(domain)):
        domain         = np.linspace(0,2*np.pi,num=10000)
        
    xa        = xc + a*np.cos(domain)
    ya        = yc + b*np.sin(domain)
    x  =  (xa)*np.cos(rot) + (ya)*np.sin(rot)
    y  = -(xa)*np.sin(rot) + (ya)*np.cos(rot) 
    return x,y
    
class RingSpreadFit:
    def __init__(self,load_step,img_num,a,b,xc,zc,
                 f,theta,coefs,intercept,fit,
                 l1_ratio,means,variances,
                 n_iter,fit_error,rel_fit_error):
                         
        self.load_step = load_step
        self.img_num = img_num        
        
        # Ellipse params
        self.a = a
        self.b = b
        self.xc = xc
        self.zc = zc
        
        # Interpolated data
        self.f = f
        self.theta = theta
        
        # Lasso fit to interpolated data
        self.coefs = coefs
        self.intercept = intercept
        self.fit = fit
        self.l1_ratio = l1_ratio
        self.means = means
        self.variances = variances
        self.n_iter = n_iter
        self.fit_error = fit_error
        self.rel_fit_error = rel_fit_error

    def plot_ellipse(self,img):
        plt.imshow(img,cmap='jet',vmin=0,vmax=200)        
        xx,zz = gen_ellipse(self.a,self.b,self.xc,self.zc,0)
        plt.plot(xx,zz,'w-')
        plt.xlim([np.min(xx),np.max(xx)])
        plt.ylim([np.min(zz),np.max(zz)])
        plt.axes().set_aspect('equal')
    
    def plot_fit(self,fig_num=0):
        plt.figure(fig_num)
        plt.plot(self.theta,self.fit,'-xr',label='Fit')
        plt.plot(self.theta,self.f,'-ob',label='Data')
        plt.ylabel('Intensity')
        plt.xlabel('Angle (Radians)')
        plt.legend()
        
    def plot_coefs(self,fig_num=0):
        plt.figure(fig_num)
        plt.plot(self.coefs,'o')
        plt.ylabel('Coefficient value')
        plt.xlabel('Coefficient index')
    
    def plot_bases_wrap(self,dtheta,fig_num=0):
        plt.figure(fig_num)        
        plt.ylabel('Intensity')
        plt.xlabel('Angle (Radians)')
         
        # Construct B
        dom = np.arange(len(self.coefs))
        keep = dom[np.abs(self.coefs)>0]
        full_fit = np.zeros(len(self.theta))

        for i in keep:
            basis = EM.gaussian_basis_wrap(len(self.f),dtheta,
                                      self.means[i],self.variances[i])
            plt.figure(fig_num)
            plt.plot(self.theta,self.coefs[i]*basis,'r-') 
        
    def plot_bases_grouped(self,dtheta,r_avg,fig_num=0):
        plt.figure(fig_num)        
        plt.ylabel('Intensity')
        plt.xlabel('Angle (Radians)')
         
        # Construct B
        dom = np.arange(len(self.coefs))
        keep = dom[np.abs(self.coefs)>0]
        last_mean = 0
        basis = np.zeros(len(self.f))
        for i in keep:
            if(self.means[i] == last_mean):
                basis += self.coefs[i]*EM.gaussian_basis_wrap(len(self.f),
                                           dtheta,
                                           self.means[i],
                                           self.variances[i])
            else:
                # Plot last basis
                plt.figure(fig_num)
                plt.plot(self.theta,basis+self.intercept)
                
                # Start a new one
                basis = self.coefs[i]*EM.gaussian_basis(len(self.f),
                                          dtheta,
                                          self.means[i],
                                          self.variances[i])
                last_mean = self.means[i]
            
              

    def scatter_amp_var(self,fig_num=0):
        plt.figure(fig_num)
        ind = self.coefs > 0
        plt.plot(self.coefs[ind],self.variances[ind],'o')
        plt.ylabel('Variance')
        plt.xlabel('Amplitude')
        
    def compute_convex_var(self):
        ind = self.coefs > 0
        self.convex_var = np.sum(self.coefs[ind]* \
                          self.variances[ind])  / \
                          np.sum(self.coefs[ind])
        
    def compute_az_var(self):
        self.az_var = np.var(self.f/np.sum(self.f))
        
    def print_params(self):
        print('Number nonzero coefficients: ' + str(np.sum(self.coefs>0)))
        print('Number of iterations: ' + str(self.n_iter))
        print('Fit error: ' + str(self.fit_error))
        print('Relative Fit error: ' + str(self.rel_fit_error))
        print('Ellipse params a='+str(self.a)+', b='+str(self.b)+\
              ', xc='+str(self.xc)+', zc='+str(self.zc) )
        print('l1_ratio = ' + str(self.l1_ratio))
        print('convex_var = ' + str(self.convex_var))
        print('az_var = ' + str(self.az_var))




def lasso_ring_fit(rsf,initial_sampling,basis_path,var_domain):
    """
    inputs:
    rsf                 RingSpreadFit object with ring image data is already loaded
    initial_sampling    integer that divides number of variances (every nth variance)
    basis_path          base file path containing basis vector files
    
    outputs:
    rsf                 RingSpreadFit object with fitted data
    
    """
    absolute_start = time.time()    
    
    print('Constructing B...')
    start = time.time()
    local_maxima = argrelmax(rsf.f)[0]
 
    
    B = np.empty([len(rsf.f),0])
    B_var = np.empty(0)
    B_mean = np.empty(0)
    for i in local_maxima:
        B_load = np.load(basis_path + str(i) + '.npy')
        B = np.hstack([B,B_load[:,::initial_sampling]])
        B_mean = np.append(B_mean, i*np.ones(len(var_domain)/initial_sampling))
        B_var  = np.append(B_var,  var_domain[::initial_sampling])
    
    timeB1 = time.time() - start
    print(timeB1)
    
    print('Lasso fitting...')
    start = time.time()
    l1_ratio = 0.08
    lasso = Lasso(alpha=l1_ratio,
                  max_iter=5000,
                  fit_intercept=0,
                  positive=True)
    # fit
    lasso.fit(B,rsf.f)        
    fit = np.dot(B,lasso.coef_) + lasso.intercept_
    fit_error     = np.linalg.norm(fit-rsf.f)
    rel_fit_error = np.linalg.norm(fit-rsf.f)/ \
                    np.linalg.norm(rsf.f)  
                    
    timeFit1 = time.time()-start
    print(timeFit1)
    
    rsf.coefs = lasso.coef_
    rsf.intercept = lasso.intercept_
    rsf.fit = fit
    rsf.l1_ratio = l1_ratio,
    rsf.means = B_mean
    rsf.variances = B_var
    rsf.n_iter = lasso.n_iter_
    rsf.fit_error = fit_error
    rsf.rel_fit_error = rel_fit_error    
    
    total_time = time.time() - absolute_start
    
    print('rel_fit_error = ' + str(rsf.rel_fit_error))
    
    return rsf, [timeB1,timeFit1,total_time]
    
def refined_lasso_ring_fit(rsf,initial_sampling,refined_sampling,basis_path,var_domain):
    """
    inputs:
    rsf                 RingSpreadFit object with ring image data is already loaded
    initial_sampling    integer that divides number of variances (every nth variance)
    refined_sampling    integer that divides number of variances (every nth variance)
    basis_path          base file path containing basis vector files
    var_domain          set of variances used to define basis functions at each mean
    
    outputs:
    rsf_refine          RingSpreadFit object with initial + refined fit
    rsf                 RingSpreadFit object with initial fit only
    benchmark           Time benchmarks for each stage of algorithm:
                        [timeB1,timeFit1,timeB2,timeFit2,total_time]  
                        Construct initial dictionary time
                        Fit using initial dictionary time
                        Construct refined dictionary time
                        Fit using refined dictionary time
                        Total run time
    
    """
    absolute_start = time.time()    
    
    print('Constructing B...')
    start = time.time()
    
    # Find local maxima
    local_maxima = argrelmax(rsf.f)[0]
 
    # Create overcomplete dictionary of Gaussians at local maxima
    B = np.empty([len(rsf.f),0])
    B_var = np.empty(0)
    B_mean = np.empty(0)
    for i in local_maxima:
        # Load precomputed matrix of Gaussians at location i
        B_load = np.load(basis_path + str(i) + '.npy')
        # Append to overcomplete dictionary
        B = np.hstack([B,B_load[:,::initial_sampling]])
        # Keep track of means and variances of these Gaussians
        B_mean = np.append(B_mean, i*np.ones(len(var_domain)/initial_sampling))
        B_var  = np.append(B_var,  var_domain[::initial_sampling])
    
    timeB1 = time.time() - start
    print(timeB1)
    
    print('Lasso fitting...')
    start = time.time()
    
    # Define a Scikit-learn Lasso solver
    l1_ratio = 0.08
    lasso = Lasso(alpha=l1_ratio,
                  max_iter=5000,
                  fit_intercept=0,
                  positive=True)
    
    # Fit data using overcomplete dictionary
    lasso.fit(B,rsf.f)        
    fit = np.dot(B,lasso.coef_) + lasso.intercept_
    fit_error     = np.linalg.norm(fit-rsf.f)
    rel_fit_error = np.linalg.norm(fit-rsf.f)/ \
                    np.linalg.norm(rsf.f)  
                    
    timeFit1 = time.time()-start
    print(timeFit1)
    
    # Save results
    rsf.coefs = lasso.coef_
    rsf.intercept = lasso.intercept_
    rsf.fit = fit
    rsf.l1_ratio = l1_ratio,
    rsf.means = B_mean
    rsf.variances = B_var
    rsf.n_iter = lasso.n_iter_
    rsf.fit_error = fit_error
    rsf.rel_fit_error = rel_fit_error

    print('Refining B...')
    start = time.time()
    
    # Compute local fit error
    loc_error = np.zeros(len(local_maxima))
    for i, max_i in enumerate(local_maxima[1:-1]):
        # Define neighborhoods as being halfway between local maxima
        up  = np.int(np.round(max_i + np.abs(local_maxima[i+1] - max_i)/2))
        dwn = np.int(np.round(max_i - np.abs(local_maxima[i-1] - max_i)/2))

        # Wrap neighborhood error computation at 0/2pi
        if(dwn > 0):
            loc_error[i] = norm(rsf.f[dwn:up]-rsf.fit[dwn:up])/norm(rsf.f)
                           
        else:
            func = np.append(rsf.f[dwn:-1],rsf.f[0:up])
            func_fit = np.append(rsf.fit[dwn:-1],rsf.fit[0:up])
            loc_error[i] = norm(func-func_fit)/norm(rsf.f)
    
    mean_error = np.mean(loc_error)
    
    # Find local maxima where error is above average 
    new_coefs = lasso.coef_
    for i, max_i in enumerate(local_maxima[loc_error > mean_error]):
        # Multiply these local maxima by 0 to remove contribution
        keep_coefs = B_mean != max_i 
        new_coefs = new_coefs*keep_coefs
        
    # Compute fit contribution from this initial stage 
    sub_fit = np.dot(B,new_coefs) + lasso.intercept_
    
    # Construct new overcomplete dictionary with only high error maxima, 
    # but with higher variance sampling density
    B_new = np.empty([len(rsf.f),0])
    B_var_new = np.empty(0)
    B_mean_new = np.empty(0)
    # Only place Gaussians where local error was above average in initial fit
    for i in local_maxima[loc_error > mean_error]:
        B_load_new = np.load(basis_path + str(i) + '.npy')
        B_new = np.hstack([B_new,B_load_new[:,::refined_sampling]])
        B_mean_new = np.append(B_mean_new, i*np.ones(len(var_domain)/refined_sampling))
        B_var_new  = np.append(B_var_new,  var_domain[::refined_sampling])

        
    timeB2 = time.time()-start
    print(timeB2)
        
    
    start = time.time()
    print('Lasso fitting refine...')
    # Fit the data only at locations where local error was high initially 
    # using identical Lasso parameters
    lasso.fit(B_new,rsf.f-sub_fit) 
    
    # Compute final fit from the intial stage and refinement stage       
    refine_fit = np.dot(B_new,lasso.coef_) + lasso.intercept_ + sub_fit
    
    # Compute fit errors
    refine_fit_error     = np.linalg.norm(refine_fit-rsf.f)
    adapt_rel_fit_error = np.linalg.norm(refine_fit-rsf.f)/ \
                          np.linalg.norm(rsf.f) 
    timeFit2 = time.time()-start
    print(timeFit2)
    
    # Store results in data structure
    rsf_refine = RingSpreadFit(rsf.load_step,rsf.img_num,
                               rsf.a,rsf.b,rsf.xc,rsf.zc,
                               rsf.f,rsf.theta,
                               0,0,0,0,0,0,0,0,0)         
    rsf_refine.coefs = np.append(new_coefs,lasso.coef_)
    rsf_refine.intercept = lasso.intercept_
    rsf_refine.fit = refine_fit
    rsf_refine.l1_ratio = l1_ratio,
    rsf_refine.means = np.append(B_mean, B_mean_new)
    rsf_refine.variances = np.append(B_var, B_var_new)
    rsf_refine.n_iter = lasso.n_iter_
    rsf_refine.fit_error = refine_fit_error
    rsf_refine.rel_fit_error = adapt_rel_fit_error
    
    total_time = time.time() - absolute_start

    benchmark = [timeB1,timeFit1,timeB2,timeFit2,total_time]    
    
    print('rel_fit_error = ' + str(rsf.rel_fit_error))
    print('rel_fit_error after adapt = ' + str(rsf_refine.rel_fit_error))
    
    return rsf_refine, rsf, benchmark
    
def generate_gaussian_basis_matrix(path,dtheta,theta0,theta1,n_var):
    
    num_theta = (theta1-theta0)/dtheta
    
    if not os.path.exists(path):
        os.mkdir(path)
        
    variances = np.linspace((dtheta),(np.pi/8),dtheta)**2
    
    B = np.empty([num_theta,0])
    
    # Generate first set of basis functions centered at 0
    for i, var in enumerate(variances):  
        basis_vec = EM.gaussian_basis_wrap(num_theta,dtheta,0,var)
        B = np.hstack([B,basis_vec[:,None]])
    
    np.save(os.path.join(path,'gauss_basis_shift_0.npy'),B)
        
    # Shift first set of basis function to create remaining sets
    for shift in range(1,num_theta):
        out_file_path = os.path.join(path,'gaus_basis_shift_' + str(shift) + '.npy')
        if not os.path.exists(out_file_path):
            print(shift)
            B_shift = np.vstack([B[-shift:,:], B[0:-shift,:]])
            np.save(out_file_path, B_shift)
            
def refined_lasso_ring_fit(rsf,initial_sampling,refined_sampling,B,var_domain):
    """
    inputs:
    rsf                 RingSpreadFit object with ring image data is already loaded
    initial_sampling    integer that divides number of variances (every nth variance)
    refined_sampling    integer that divides number of variances (every nth variance)
    basis_path          base file path containing basis vector files
    var_domain          set of variances used to define basis functions at each mean
    
    outputs:
    rsf_refine          RingSpreadFit object with initial + refined fit
    rsf                 RingSpreadFit object with initial fit only
    benchmark           Time benchmarks for each stage of algorithm:
                        [timeB1,timeFit1,timeB2,timeFit2,total_time]  
                        Construct initial dictionary time
                        Fit using initial dictionary time
                        Construct refined dictionary time
                        Fit using refined dictionary time
                        Total run time
    
    """
    absolute_start = time.time()    
    
    print('Constructing B...')
    start = time.time()
    
    # Find local maxima
    local_maxima = argrelmax(rsf.f)[0]
    
    print('Lasso fitting...')
    start = time.time()
    
    # Define a Scikit-learn Lasso solver
    l1_ratio = 0.08
    lasso = Lasso(alpha=l1_ratio,
                  max_iter=5000,
                  fit_intercept=0,
                  positive=True)
    
    # Fit data using overcomplete dictionary
    lasso.fit(B,rsf.f)        
    fit = np.dot(B,lasso.coef_) + lasso.intercept_
    fit_error     = np.linalg.norm(fit-rsf.f)
    rel_fit_error = np.linalg.norm(fit-rsf.f)/ \
                    np.linalg.norm(rsf.f)  
                    
    timeFit1 = time.time()-start
    print(timeFit1)
    
    # Save results
    rsf.coefs = lasso.coef_
    rsf.intercept = lasso.intercept_
    rsf.fit = fit
    rsf.l1_ratio = l1_ratio,
    rsf.means = B_mean
    rsf.variances = B_var
    rsf.n_iter = lasso.n_iter_
    rsf.fit_error = fit_error
    rsf.rel_fit_error = rel_fit_error

    print('Refining B...')
    start = time.time()
    
    # Compute local fit error
    loc_error = np.zeros(len(local_maxima))
    for i, max_i in enumerate(local_maxima[1:-1]):
        # Define neighborhoods as being halfway between local maxima
        up  = np.int(np.round(max_i + np.abs(local_maxima[i+1] - max_i)/2))
        dwn = np.int(np.round(max_i - np.abs(local_maxima[i-1] - max_i)/2))

        # Wrap neighborhood error computation at 0/2pi
        if(dwn > 0):
            loc_error[i] = norm(rsf.f[dwn:up]-rsf.fit[dwn:up])/norm(rsf.f)
                           
        else:
            func = np.append(rsf.f[dwn:-1],rsf.f[0:up])
            func_fit = np.append(rsf.fit[dwn:-1],rsf.fit[0:up])
            loc_error[i] = norm(func-func_fit)/norm(rsf.f)
    
    mean_error = np.mean(loc_error)
    
    # Find local maxima where error is above average 
    new_coefs = lasso.coef_
    for i, max_i in enumerate(local_maxima[loc_error > mean_error]):
        # Multiply these local maxima by 0 to remove contribution
        keep_coefs = B_mean != max_i 
        new_coefs = new_coefs*keep_coefs
        
    # Compute fit contribution from this initial stage 
    sub_fit = np.dot(B,new_coefs) + lasso.intercept_
    
    # Construct new overcomplete dictionary with only high error maxima, 
    # but with higher variance sampling density
    B_new = np.empty([len(rsf.f),0])
    B_var_new = np.empty(0)
    B_mean_new = np.empty(0)
    # Only place Gaussians where local error was above average in initial fit
    for i in local_maxima[loc_error > mean_error]:
        B_load_new = np.load(basis_path + str(i) + '.npy')
        B_new = np.hstack([B_new,B_load_new[:,::refined_sampling]])
        B_mean_new = np.append(B_mean_new, i*np.ones(len(var_domain)/refined_sampling))
        B_var_new  = np.append(B_var_new,  var_domain[::refined_sampling])

        
    timeB2 = time.time()-start
    print(timeB2)
        
    
    start = time.time()
    print('Lasso fitting refine...')
    # Fit the data only at locations where local error was high initially 
    # using identical Lasso parameters
    lasso.fit(B_new,rsf.f-sub_fit) 
    
    # Compute final fit from the intial stage and refinement stage       
    refine_fit = np.dot(B_new,lasso.coef_) + lasso.intercept_ + sub_fit
    
    # Compute fit errors
    refine_fit_error     = np.linalg.norm(refine_fit-rsf.f)
    adapt_rel_fit_error = np.linalg.norm(refine_fit-rsf.f)/ \
                          np.linalg.norm(rsf.f) 
    timeFit2 = time.time()-start
    print(timeFit2)
    
    # Store results in data structure
    rsf_refine = RingSpreadFit(rsf.load_step,rsf.img_num,
                               rsf.a,rsf.b,rsf.xc,rsf.zc,
                               rsf.f,rsf.theta,
                               0,0,0,0,0,0,0,0,0)         
    rsf_refine.coefs = np.append(new_coefs,lasso.coef_)
    rsf_refine.intercept = lasso.intercept_
    rsf_refine.fit = refine_fit
    rsf_refine.l1_ratio = l1_ratio,
    rsf_refine.means = np.append(B_mean, B_mean_new)
    rsf_refine.variances = np.append(B_var, B_var_new)
    rsf_refine.n_iter = lasso.n_iter_
    rsf_refine.fit_error = refine_fit_error
    rsf_refine.rel_fit_error = adapt_rel_fit_error
    
    total_time = time.time() - absolute_start

    benchmark = [timeB1,timeFit1,timeB2,timeFit2,total_time]    
    
    print('rel_fit_error = ' + str(rsf.rel_fit_error))
    print('rel_fit_error after adapt = ' + str(rsf_refine.rel_fit_error))
    
    return rsf_refine, rsf, benchmark
    
class Ring2DFit:
    def __init__(self,load_step,img_num,a,b,xc,yc,
                 f,mask_xy,img_shape,
                 az_means,az_variances,rad_dx,rad_dy,rad_variances,
                 coefs,fit,l1_ratio,n_iter,fit_error,rel_fit_error):
                         
        self.load_step = load_step
        self.img_num = img_num        
        
        # Ellipse params
        self.a = a
        self.b = b
        self.xc = xc
        self.yc = yc
        
        # Support of ring image
        self.f = f
        self.mask_xy = mask_xy
        self.img_shape = img_shape
        
        # Lasso fit to interpolated data
        self.coefs = coefs

        self.fit = fit
        self.l1_ratio = l1_ratio
        self.az_means = az_means
        self.az_variances = az_variances
        self.rad_dx = rad_dx
        self.rad_dy = rad_dy
        self.rad_variances = rad_variances
        self.n_iter = n_iter
        self.fit_error = fit_error
        self.rel_fit_error = rel_fit_error

    def plot_ellipse(self,img):
        plt.imshow(img,cmap='jet',vmin=0,vmax=200)        
        xx,zz = gen_ellipse(self.a,self.b,self.xc,self.zc,0)
        plt.plot(xx,zz,'w-')
        plt.xlim([np.min(xx),np.max(xx)])
        plt.ylim([np.min(zz),np.max(zz)])
        plt.axes().set_aspect('equal')
    
    def display_fit_image(self,fig_num=0):
        fit_img = np.zeros(self.img_shape)
        fit_img[self.mask_xy[1],self.mask_xy[0]] = self.fit        
        plt.figure(fig_num)
        plt.imshow(fit_img,cmap='jet',vmin=0,vmax=200,interpolation='nearest')   
        plt.colorbar()  
        
    def display_f_image(self,fig_num=0):
        fit_img = np.zeros(self.img_shape)
        fit_img[self.mask_xy[1],self.mask_xy[0]] = self.f        
        plt.figure(fig_num)
        plt.imshow(fit_img,cmap='jet',vmin=0,vmax=200,interpolation='nearest')   
        plt.colorbar()
    
    def print_params(self):
        print('Number nonzero coefficients: ' + str(np.sum(self.coefs>0)))
        print('Number of iterations: ' + str(self.n_iter))
        print('Fit error: ' + str(self.fit_error))
        print('Relative Fit error: ' + str(self.rel_fit_error))
        print('Ellipse params a='+str(self.a)+', b='+str(self.b)+\
              ', xc='+str(self.xc)+', zc='+str(self.zc) )
        print('l1_ratio = ' + str(self.l1_ratio))
        print('convex_var = ' + str(self.convex_var))
        print('az_var = ' + str(self.az_var))
        
def lasso_ring_2Dfit(r2Df,B):
    """
    inputs:
    rsf                 RingSpreadFit object with ring image data is already loaded
    initial_sampling    integer that divides number of variances (every nth variance)
    refined_sampling    integer that divides number of variances (every nth variance)
    basis_path          base file path containing basis vector files
    var_domain          set of variances used to define basis functions at each mean
    
    outputs:
    rsf_refine          RingSpreadFit object with initial + refined fit
    rsf                 RingSpreadFit object with initial fit only
    benchmark           Time benchmarks for each stage of algorithm:
                        [timeB1,timeFit1,timeB2,timeFit2,total_time]  
                        Construct initial dictionary time
                        Fit using initial dictionary time
                        Construct refined dictionary time
                        Fit using refined dictionary time
                        Total run time
    
    """
    absolute_start = time.time()    
    
    print('Lasso fitting...')
    start = time.time()
    
    # Define a Scikit-learn Lasso solver
    l1_ratio = 0.08
    lasso = Lasso(alpha=l1_ratio,
                  max_iter=5000,
                  fit_intercept=0,
                  positive=True)
    
    # Fit data using overcomplete dictionary
    lasso.fit(B,r2Df.f)        
    fit = np.dot(B,lasso.coef_)
    fit_error     = np.linalg.norm(fit-r2Df.f)
    rel_fit_error = np.linalg.norm(fit-r2Df.f)/ \
                    np.linalg.norm(r2Df.f)  
                    
    timeFit1 = time.time()-start
    print(timeFit1)
    
    # Save results
    r2Df.coefs = lasso.coef_
    r2Df.fit = fit
    r2Df.l1_ratio = l1_ratio,
    r2Df.n_iter = lasso.n_iter_
    r2Df.fit_error = fit_error
    r2Df.rel_fit_error = rel_fit_error
    
    benchmark = timeFit1    
    
    print('rel_fit_error = ' + str(r2Df.rel_fit_error))
    
    return r2Df, benchmark

def circulant_lasso_ring_2Dfit(r2Df,B,B0,B0_tran):
    """
    inputs:
    rsf                 RingSpreadFit object with ring image data is already loaded
    initial_sampling    integer that divides number of variances (every nth variance)
    refined_sampling    integer that divides number of variances (every nth variance)
    basis_path          base file path containing basis vector files
    var_domain          set of variances used to define basis functions at each mean
    
    outputs:
    rsf_refine          RingSpreadFit object with initial + refined fit
    rsf                 RingSpreadFit object with initial fit only
    benchmark           Time benchmarks for each stage of algorithm:
                        [timeB1,timeFit1,timeB2,timeFit2,total_time]  
                        Construct initial dictionary time
                        Fit using initial dictionary time
                        Construct refined dictionary time
                        Fit using refined dictionary time
                        Total run time
    
    """
    absolute_start = time.time()    
    
    print('Lasso fitting...')
    start = time.time()
    
    # Define a Scikit-learn Lasso solver
    l1_ratio = 0.08
    max_iters = 5000
    lasso = Lasso(alpha=l1_ratio,
                  fit_intercept=0,
                  positive=True)
    
    # Fit data using overcomplete dictionary
    coefs = Lasso.circulant_fista.fit(B,B0,B0_tran,maxima,r2Df.f,  
                                      l1_ratio, max_iters, 
                                      positive=1, benchmark=0)    
    fit = np.dot(B,lasso.coef_)
    fit_error     = np.linalg.norm(fit-r2Df.f)
    rel_fit_error = np.linalg.norm(fit-r2Df.f)/ \
                    np.linalg.norm(r2Df.f)  
                    
    timeFit1 = time.time()-start
    print(timeFit1)
    
    # Save results
    r2Df.coefs = lasso.coef_
    r2Df.fit = fit
    r2Df.l1_ratio = l1_ratio,
    r2Df.n_iter = lasso.n_iter_
    r2Df.fit_error = fit_error
    r2Df.rel_fit_error = rel_fit_error
    
    benchmark = timeFit1    
    
    print('rel_fit_error = ' + str(r2Df.rel_fit_error))
    
    return r2Df, benchmark


def soft_thresh(rho, lam, z):
    if(rho < -lam/2):
        y = (rho + lam/2)/z
    elif(rho > lam/2):
        y = (rho - lam/2)/z  
    else:
        y = 0
    return y
        

def lasso_solver(A,b,lam,max_iter):
    m,n = A.shape
    
    # precompute sum of sqaures of columns
    z = np.sum(A**2,0)    
    
    # initialize x
    x = np.zeros((n,1))

    not_converged = 1
    iters = 0
    while(not_converged):
        
        print(iters)
        x_old = x.copy()
        for j in range(n):
            rho = np.inner(A[:,j], b-A.dot(x).ravel())
        for j in range(n):
            x[j] = 0
            x[j] = soft_thresh(rho, lam, z[j])
        
        iters += 1
        print(norm(x - x_old))
        not_converged = norm(x - x_old) > 10**(-8)
        not_converged = iters < max_iter
            
    fit = A.dot(x)
    rel_error = norm( b - fit )/norm(b)
    return x, fit, rel_error
    