# -*- coding: utf-8 -*-
"""
L-BFGS-B Quick Test Code

Created on Tue Feb 21 10:28:28 2017

@author: dbanco02
"""

from __future__ import print_function
from __future__ import division

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from numpy.linalg import norm 

A = np.random.rand(200,4000)
b = np.random.rand(200,)
x0 = np.ones((4000,))

def objective_func(x,A,b):
    objective = np.sum((A.dot(x)-b)**2) + np.sum(np.abs(x))
    return objective
 
def gradient_func(x,A,b):
    gradient = 2*A.T.dot(A.dot(x)-b) + 2*x/np.sqrt(x**2 + 10**(-8))
    return gradient

callbackf = lambda x: print(np.sum((A.dot(x)-b)**2) + np.sum(np.abs(x)))
f_func = lambda x: np.sum((A.dot(x)-b)**2) + np.sum(np.abs(x))
g_func = lambda x: 2*A.T.dot(A.dot(x)-b) + 2*x/np.sqrt(x**2 + 10**(-8))

x_bar1 = opt.fmin_l_bfgs_b(func=f_func,
                          x0=x0,
                          fprime = g_func,
                          iprint=1,
                          callback=callbackf,
                          disp=1)

x_bar2 = opt.fmin_l_bfgs_b(func=objective_func,
                          x0=x0,
                          fprime = gradient_func,
                          args=(A,b),
                          iprint=1,
                          disp=1)

f1 = A.dot(x_bar1[0])
f2 = A.dot(x_bar2[0])

print(norm(f1-b)/norm(b))
print(norm(f2-b)/norm(b))
plt.plot(b)
plt.plot(f1)
plt.plot(f2)


