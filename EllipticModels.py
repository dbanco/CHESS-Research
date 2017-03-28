# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 18:19:23 2016

Ellipse models


@author: Dan
"""
import numpy as np

"""
Generates an elliptical basis defined on a grid. This is most useful for 
viewing the generated basis functions

EX:

## Points where we have data for a single ring
x = np.linspace(-1.05,1.05,200)
z = np.linspace(-1.05,1.05,200)

# Eccentricity of circle fit
r = 1
delta = 0
a = r - delta
b = r + delta

# Rotation of ellipse
theta = np.pi/4

# Number of points sampled along ellipse
numt = 360

# Place Gaussian every tskip points
tskip = 10
numGauss = int(np.floor(numt/tskip))

# Variance perpendicular to the curve
# fit this separately
ssq_perp = (5e-3)**2 
   
# Variance parallel to the curve 
# Should not really go below perpendicular variance
# Use multi-scale processing theory to select variances
   
ssq_prll = (8e-2)**2

B, x0, y0 = ellipse_basis_grid(x,z,a,b,theta,numt,tskip,ssq_perp,ssq_prll[i])

plt.figure(1)
plt.contourf(x,z,np.sum(B,axis=2))
plt.colorbar()

"""

def ellipse_basis_grid(x,y,a,b,theta,numt,tskip,ssq_perp,ssq_prll):
    """
    inputs:
        x,y: grid of points over which basis functions will be evaluated
        a: half x-axis length of ellipse
        b: half y-axis length of ellipse
        theta: tilt of ellipse in radians
        numt: number of points to be sampled around the ellipse
        tskip: will create basis functions every tskip points in t
        ssq_perp: Variance perpendicular to the curve
        ssq_prll: Variance paralell to the curve    
        
    outputs:
        B: 3D set of basis functions
        xe,ye: coordinates of ellipse      """
    
    
    ct = np.cos(theta)
    st = np.sin(theta)
    
    # Sample the ellipse at numt points
    t = np.linspace(0,1,numt)
    dt = t[1] - t[0]
    
    # Build the x and y coordinates of the ellipse
    x0 = a*ct*np.cos(2*np.pi*t) - b*st*np.sin(2*np.pi*t)
    y0 = a*st*np.cos(2*np.pi*t) + b*ct*np.sin(2*np.pi*t)
    x0 = x0.ravel()
    y0 = y0.ravel()
    
    # Find the length of the ellipse.  Probably there is a formula for this,
    # but if we want to change the curve, then will need this code anyway.
    # Basic idea: approximate legnth as the sum of the lengths between samples
    # along the curve.
    curvelength = 0
    for idx in range(1,numt):
        curvelength += np.sqrt( (x0[idx]-x0[idx-1])**2 +
                                (y0[idx]-y0[idx-1])**2 )
                                
    # Multiply by dt to make this a real length (i.e., approximation the true
    # value of the continuous integral.
    curvelength = curvelength*dt
    
    # Indices where we will generate basis functions
    tidx_all = range(0,numt-1,tskip)
    
    # Hold the "image" of our basis function
    B = np.zeros([len(x),len(y),len(tidx_all)])
 
    for tidx in range(len(tidx_all)):
        tidx0 = tidx_all[tidx]
        print('point ' + str(tidx) + ' of ' + str(len(tidx_all)))
        # Loop over all points in grid
        for xidx in range(len(x)):
            for yidx in range(len(y)):
                # For any x and y, find the distance to the curve and the 
                # t index of the point on the curve for which the distance is 
                # minimized. Note that ths distance we find will be the 
                # perpendicular distance to the curve.
                dist = (x[xidx]-x0)**2 + (y[yidx]-y0)**2
                dsq_perp = np.amin(dist)
                tidx1    = np.argmin(dist)
                
                # Knowing the point on the curve where the basis function is 
                # to be centered and the point on the curve closest to (x,y), 
                # we compute the distance between the two.
            
                # Start at the point with the smallest index and end at the 
                # point with the largest                
                tidxa = np.amin([tidx0,tidx1]) 
                tidxb = np.amax([tidx0,tidx1])
                
                # Again, add up the distances between all points between the two...
                d_prll = 0
                for idx in range(tidxa,tidxb+1):
                    
                    d_prll += np.sqrt((x0[idx]-x0[idx-1])**2 + (y0[idx]-y0[idx-1])**2);
                
                # ...and scale by dt
                d_prll = d_prll*dt
                
                # Since the curve is closed, there are two ways we can go from one
                # point to the other.  We want the shorter of the two.  If the
                # distance we computed was more than half the curvlength, we must
                # have gone around the wrong way.  This is easily fixed though
                # since the right way will be the length of the curve minus the
                # length of our wrong path.
                if (d_prll > curvelength/2):
                    d_prll = curvelength - d_prll
                
                dsq_prll = d_prll**2
                
                # So now we have the parallel and perpendicular distances.  We put
                # them together to build the basis vector.
                B[xidx,yidx,tidx] = np.exp(-(dsq_perp/ssq_perp + dsq_prll/ssq_prll))
                
    return B, x0, y0

"""
Generates an ellipse basis for a set of coordinates pos of the form (x,y)
where x and y are numpy arrays of shape (n,). This is useful for fitting basis
functions to particular set of data points

EX:

See example for ellipse_basis_grid above. Define pos = (x,y). Note that this
will only evaluate the basis at these exact x,y coordinates instead of using 
the x,y coordinates to fill a grid.
"""

def ellipse_basis(pos,a,b,theta,numt,tskip,ssq_perp,ssq_prll):
    """
    inputs:
        x,y: grid of points over which basis functions will be evaluated
        a: half x-axis length of ellipse
        b: half y-axis length of ellipse
        theta: tilt of ellipse in radians
        numt: number of points to be sampled around the ellipse
        tskip: will create basis functions every tskip points in t
        ssq_perp: Variance perpendicular to the curve
        ssq_prll: Variance paralell to the curve    
        
    outputs:
        B: 3D set of basis functions
        xe,ye: coordinates of ellipse      """
    
    
    ct = np.cos(theta)
    st = np.sin(theta)
    
    # Sample the ellipse at numt points
    t = np.linspace(0,1,numt)
    dt = t[1] - t[0]
    
    # Build the x and y coordinates of the ellipse
    x0 = a*ct*np.cos(2*np.pi*t) - b*st*np.sin(2*np.pi*t)
    y0 = a*st*np.cos(2*np.pi*t) + b*ct*np.sin(2*np.pi*t)
    x0 = x0.ravel()
    y0 = y0.ravel()
    
    # Find the length of the ellipse.  Probably there is a formula for this,
    # but if we want to change the curve, then will need this code anyway.
    # Basic idea: approximate legnth as the sum of the lengths between samples
    # along the curve.
    curvelength = 0
    for idx in range(1,numt):
        curvelength += np.sqrt( (x0[idx]-x0[idx-1])**2 +
                                (y0[idx]-y0[idx-1])**2 )
                                
    # Multiply by dt to make this a real length (i.e., approximation the true
    # value of the continuous integral.
    curvelength = curvelength*dt
    
    # Indices where we will generate basis functions
    tidx_all = range(0,numt-1,tskip)
    
    # Hold the "image" of our basis function
    x = pos[0]
    y = pos[1]
    B = np.zeros([len(x),len(tidx_all)])
 
    for tidx in range(len(tidx_all)):
        tidx0 = tidx_all[tidx]
        print('point ' + str(tidx) + ' of ' + str(len(tidx_all)))
        # Loop over all points in grid
        for posidx in range(len(x)):
            # For any x and y, find the distance to the curve and the 
            # t index of the point on the curve for which the distance is 
            # minimized. Note that ths distance we find will be the 
            # perpendicular distance to the curve.
            dist = (x[posidx]-x0)**2 + (y[posidx]-y0)**2
            dsq_perp = np.amin(dist)
            tidx1    = np.argmin(dist)
            
            # Knowing the point on the curve where the basis function is 
            # to be centered and the point on the curve closest to (x,y), 
            # we compute the distance between the two.
        
            # Start at the point with the smallest index and end at the 
            # point with the largest                
            tidxa = np.amin([tidx0,tidx1]) 
            tidxb = np.amax([tidx0,tidx1])
            
            # Again, add up the distances between all points between the two...
            d_prll = 0
            for idx in range(tidxa,tidxb+1):
                
                d_prll += np.sqrt((x0[idx]-x0[idx-1])**2 + (y0[idx]-y0[idx-1])**2);
            
            # ...and scale by dt
            d_prll = d_prll*dt
            
            # Since the curve is closed, there are two ways we can go from one
            # point to the other.  We want the shorter of the two.  If the
            # distance we computed was more than half the curvlength, we must
            # have gone around the wrong way.  This is easily fixed though
            # since the right way will be the length of the curve minus the
            # length of our wrong path.
            if (d_prll > curvelength/2):
                d_prll = curvelength - d_prll
            
            dsq_prll = d_prll**2
            
            # So now we have the parallel and perpendicular distances.  We put
            # them together to build the basis vector.
            B[posidx,tidx] = np.exp(-(dsq_perp/ssq_perp + dsq_prll/ssq_prll))
            
    return B, x0, y0

def ellipse_basis_single_spot(pos,a,b,xc,yc,dx,dy,theta,num_theta,var_perp,var_prll):
    """
    inputs:
        mask: points over which basis functions will be evaluated
        pos: points for which elliptic gaussian is defined
        a: half x-axis length of ellipse
        b: half y-axis length of ellipse
        xc: x-axis center of ellipse
        yc: y-axis center of ellipse
        dx: x-axis offset
        dy: y-axis offset
        theta: azimuthal mean of spot
        num_theta: number of azimuthal samples
        ssq_perp: Variance perpendicular to the curve
        ssq_prll: Variance paralell to the curve    
        
    outputs:
        B: 3D set of basis functions
        xe,ye: coordinates of ellipse      """
    
    # Sample the ellipse at numt points
    t = np.linspace(0,2*np.pi,num_theta)
    dtheta = t[1] - t[0]
    
    # Build the x and y coordinates of the ellipse
    x0 = a*np.cos(t) + xc - dx
    y0 = b*np.sin(t) + yc - dy
    x0 = x0.ravel()
    y0 = y0.ravel()    
    
    # Hold the "image" of our basis function
    x = pos[0]
    y = pos[1]    
    

    # Find the length of the ellipse.  Probably there is a formula for this,
    # but if we want to change the curve, then will need this code anyway.
    # Basic idea: approximate legnth as the sum of the lengths between samples
    # along the curve.
    curvelength = 2*np.pi
    
    tidx0 = np.argmin(abs(t - theta))
    B = np.zeros([len(x)])
 
    # Loop over all points in grid
    for posidx in range(len(x)):
        # For any x and y, find the distance to the curve and the 
        # t index of the point on the curve for which the distance is 
        # minimized. Note that ths distance we find will be the 
        # perpendicular distance to the curve.
        dist = (x[posidx]-x0)**2 + (y[posidx]-y0)**2
        dsq_perp = np.amin(dist)
        tidx1    = np.argmin(dist)
        
        # Knowing the point on the curve where the basis function is 
        # to be centered and the point on the curve closest to (x,y), 
        # we compute the distance between the two.
    
        # Start at the point with the smallest index and end at the 
        # point with the largest                
        tidxa = np.amin([tidx0,tidx1]) 
        tidxb = np.amax([tidx0,tidx1])
        
        # Again, add up the distances between all points between the two...
       
        d_prll = (tidxb-tidxa)*dtheta
        
        # Since the curve is closed, there are two ways we can go from one
        # point to the other.  We want the shorter of the two.  If the
        # distance we computed was more than half the curvlength, we must
        # have gone around the wrong way.  This is easily fixed though
        # since the right way will be the length of the curve minus the
        # length of our wrong path.
        if (d_prll > curvelength/2):
            d_prll = curvelength - d_prll
        
        dsq_prll = d_prll**2
        
        # So now we have the parallel and perpendicular distances.  We put
        # them together to build the basis vector.
        B[posidx] = np.exp(-(dsq_perp/var_perp + dsq_prll/var_prll))
                             # pixel units       # radian units
        
    return B


"""
#%% Test Gaussian on ellipse basis

# Setup basic coordinate system.  Both x and y run from -1 to 1.  Could
# easily change this as needed.  Also could make sampling more or less
# dense as needed.
x = np.linspace(-1,1,20)
y = np.linspace(-1,1,20)

# Look at an ellipse
a = 0.7 # half the length of the x axis
b = 0.75 # half the length of the y axis
theta = np.pi/10 # Tilt angle in radians

numt = 100 # number of points to be sampled around the ellipse
tskip = 10 # will create basis functions every tskip points in t

ssq_perp = (0.03)**2    # Variance perpendicular to the curve
ssq_prll = (0.04/20)**2 # Variance parallel to the curve

# Generate the basis
[B,xe,ye] = ellipse_basis(x,y,a,b,theta,numt,tskip,ssq_perp,ssq_prll);

# Plot the "sum" of all the basis elements
plt.contourf(x,y,np.sum(B,axis=2))
# Overlay the points on the elipse just to make sure we are OK
plt.plot(xe,ye,'x')
"""


def gaussian_basis(num_x,dx,r,mean,variance):
    """
    inputs:
            num_x      number of data points
            dx         distance between data points
            mean       mean of gaussian function
            variance   variance of gaussian function
       
    outputs:
            B           basis vector
    """
    
    idx = np.arange(0,num_x)
    opposite = (idx[mean-np.floor(num_x/2)] + idx[mean-np.ceil(num_x/2)])/2
    dist1 = np.abs(mean - idx)
    dist2 = num_x/2 - np.abs(opposite - idx)
    dist = np.minimum(dist1,dist2)*dx*r
    # This distance may need to be multipled by dtheta to units of theta 
    
    # Square the distance
    dist_sq = dist**2
    
    # Build basis vector
    B = np.exp(-dist_sq/variance)
            
    return B
    
def gaussian_basis_wrap(num_theta,dtheta,mean,variance):
    """
    inputs:
            num_theta      number of data points
            dtheta         distance between data points
            mean       mean of gaussian function
            variance   variance of gaussian function
       
    outputs:
            B           basis vector
    """
    
    idx = np.arange(0,num_theta)
    opposite = (idx[int(mean-np.floor(num_theta/2))] + idx[int(mean-np.ceil(num_theta/2))])/2
    dist1 = np.abs(mean - idx)
    dist2 = num_theta/2 - np.abs(opposite - idx)
    dist = np.minimum(dist1,dist2)*dtheta

    # Square the distance
    dist_sq = dist**2
    
    # Build basis vector
    B = np.exp(-dist_sq/variance)
            
    return B

def gaussian_basis_wrap_2D(num_theta,dtheta,m_theta,v_theta,num_r,dr,m_r,v_r):
    """
    inputs:
            num_theta      number of data points
            dtheta         distance between data points
            mean       mean of gaussian function
            variance   variance of gaussian function
       
    outputs:
            B           basis vector
    """
    
    # Compute theta distances with wrapping
    idx = np.arange(0,num_theta)
    opposite = (idx[int(m_theta-np.floor(num_theta/2))] + idx[int(m_theta-np.ceil(num_theta/2))])/2
    dist1 = np.abs(m_theta - idx)
    dist2 = num_theta/2 - np.abs(opposite - idx)
    dist = np.minimum(dist1,dist2)*dtheta
    dist_sq_theta = dist**2    # num_theta length vector
    
    # Compute radial distance with wrapping
    idx = np.arange(0,num_r)
    opposite = (idx[int(m_r-np.floor(num_r/2))] + idx[int(m_r-np.ceil(num_r/2))])/2
    dist1 = np.abs(m_r - idx)
    dist2 = num_r/2 - np.abs(opposite - idx)
    dist = np.minimum(dist1,dist2)*dr
    dist_sq_r = dist**2    # num_r length vector

    all_dist_sq_theta = dist_sq_theta.copy()
    all_dist_sq_r = dist_sq_r[0]*np.ones((num_theta,))
    for i in range(1,num_r):
        all_dist_sq_theta = np.append(all_dist_sq_theta,dist_sq_theta)
        all_dist_sq_r = np.append(all_dist_sq_r, dist_sq_r[i]*np.ones((num_theta,)))

    # Build basis vector
    B = np.exp(-all_dist_sq_theta/v_theta - all_dist_sq_r/v_r)
            
    return B

def gaussian_basis_wrap_2D_tube(num_theta,dtheta,v_theta,num_rad,dr,v_r):
    """
    inputs:
            num_theta      number of data points
            dtheta         distance between data points
            mean       mean of gaussian function
            variance   variance of gaussian function
       
    outputs:
            B           basis vector
    """
    
    # Compute theta distances with wrapping
    m_theta = np.arange(0,num_theta)
    opposite = m_theta + num_theta/2
    opposite[m_theta>=num_theta/2] = np.mod(m_theta[m_theta>=num_theta/2] + num_theta/2,num_theta/2)
    dist1 = np.abs(m_theta)
    dist2 = num_theta/2 - np.abs(opposite)
    dist = np.minimum(dist1,dist2)*dtheta
    dist_sq_theta = dist**2    # num_theta length vector
    
    
    # Compute radial distance with wrapping
    m_rad = np.arange(0,num_rad)
    opposite = m_rad + num_rad/2
    opposite[m_rad>=num_rad/2] = np.mod(m_rad[m_rad>=num_rad/2] + num_rad/2,num_rad/2)
    dist1 = np.abs(m_rad)
    dist2 = num_rad/2 - np.abs(opposite)
    dist = np.minimum(dist1,dist2)*dr
    dist_sq_r = dist**2    # num_r length vector

    all_dist_sq_theta = dist_sq_theta.copy()
    all_dist_sq_r = dist_sq_r[0]*np.ones((num_theta,))
    coor_theta = m_theta.copy()
    coor_rad = m_rad[0]*np.ones((num_theta,))
    for i in range(1,num_rad):
        all_dist_sq_theta = np.append(all_dist_sq_theta,dist_sq_theta)
        all_dist_sq_r = np.append(all_dist_sq_r, dist_sq_r[i]*np.ones((num_theta,)))
        coor_theta = np.append(coor_theta,m_theta)
        coor_rad = np.append(coor_rad,m_rad[i]*np.ones((num_theta,)))
    # Build basis vector
    B = np.exp(-all_dist_sq_theta/v_theta - all_dist_sq_r/v_r)
    domain = np.array((coor_theta,coor_rad))
            
    return B, domain

def gaussian_basis_wrap_2D_shift_tube(num_theta,dtheta,v_theta,num_rad,dr,v_r,idx_t,idx_r):
    """
    inputs:
            num_theta      number of data points
            dtheta         distance between data points
            mean       mean of gaussian function
            variance   variance of gaussian function
       
    outputs:
            B           basis vector
    """
    
    # Compute theta distances with wrapping
    m_theta = np.arange(0,num_theta)
    opposite = m_theta + num_theta/2
    opposite[m_theta>=num_theta/2] = np.mod(m_theta[m_theta>=num_theta/2] + num_theta/2,num_theta/2)
    dist1 = np.abs(m_theta - idx_t)
    dist2 = num_theta/2 - np.abs(opposite - idx_t)
    dist = np.minimum(dist1,dist2)*dtheta
    dist_sq_theta = dist**2    # num_theta length vector
    
    
    # Compute radial distance with wrapping
    m_rad = np.arange(0,num_rad)
    opposite = m_rad + num_rad/2
    opposite[m_rad>=num_rad/2] = np.mod(m_rad[m_rad>=num_rad/2] + num_rad/2,num_rad/2)
    dist1 = np.abs(m_rad - idx_r)
    dist2 = num_rad/2 - np.abs(opposite - idx_r)
    dist = np.minimum(dist1,dist2)*dr
    dist_sq_r = dist**2    # num_r length vector

    all_dist_sq_theta = dist_sq_theta.copy()
    all_dist_sq_r = dist_sq_r[0]*np.ones((num_theta,))
    coor_theta = m_theta.copy()
    coor_rad = m_rad[0]*np.ones((num_theta,))
    for i in range(1,num_rad):
        all_dist_sq_theta = np.append(all_dist_sq_theta,dist_sq_theta)
        all_dist_sq_r = np.append(all_dist_sq_r, dist_sq_r[i]*np.ones((num_theta,)))
        coor_theta = np.append(coor_theta,m_theta)
        coor_rad = np.append(coor_rad,m_rad[i]*np.ones((num_theta,)))
    # Build basis vector
    B = np.exp(-all_dist_sq_theta/v_theta - all_dist_sq_r/v_r)
    domain = np.array((coor_theta,coor_rad))
            
    return B, domain

def unshifted_basis_matrix_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad):
    A0_stack = np.zeros((num_rad,num_theta,var_theta.shape[0],var_rad.shape[0]))
    for t in range(var_theta.shape[0]):
        for r in range(var_rad.shape[0]):
            A0 = gaussian_basis_wrap_2D(num_theta,dtheta,  0,  var_theta[t], 
                                        num_rad,  drad,    0,  var_rad[r])
            A0_stack[:,:,t,r] = A0.reshape((num_rad,num_theta))
    return A0_stack
    
def unshifted_basis_matrix_ft_stack(var_theta,var_rad,dtheta,drad,num_theta,num_rad):
    A0ft_stack = np.zeros((num_rad,num_theta,var_theta.shape[0],var_rad.shape[0]),dtype=complex)
    for t in range(var_theta.shape[0]):
        for r in range(var_rad.shape[0]):
            A0 = gaussian_basis_wrap_2D(num_theta,dtheta,  0,  var_theta[t], 
                                        num_rad,  drad,    0,  var_rad[r])
            A0ft_stack[:,:,t,r] = np.fft.fft2(A0.reshape((num_rad,num_theta)))
    return A0ft_stack

def unshifted_basis_matrix_ft_list(var_theta,var_rad,dtheta,drad,num_theta,num_rad):
    A0ft_list = []
    for t in range(var_theta.shape[0]):
        for r in range(var_rad.shape[0]):
            A0 = gaussian_basis_wrap_2D(num_theta,dtheta,  0,  var_theta[t], 
                                        num_rad,  drad,    0,  var_rad[r])
            A0ft_list.append( [np.fft.fft2(A0.reshape((num_rad,num_theta)))] )
    return A0ft_list
            
def unshifted_basis_ft_svd_list(var_theta,var_rad,dtheta,drad,num_theta,num_rad):
    A0ft_list = []
    for t in range(var_theta.shape[0]):
        for r in range(var_rad.shape[0]):
            A0 = gaussian_basis_wrap_2D(num_theta,dtheta,  0,  var_theta[t], 
                                        num_rad,  drad,    0,  var_rad[r])
            u,s,v = np.linalg.svd(np.fft.fft2(A0.reshape((num_rad,num_theta))))
            A0ft_list.append([ u[:,0], s[0], v[0,:] ])
    return A0ft_list
        
