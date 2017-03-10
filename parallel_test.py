# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 17:26:18 2017

@author: dbanco02
"""
import numpy as np
import multiprocessing as mp
from functools import partial
import time

def f1(x,A):
    return [ np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x)),np.fft.fft(A.dot(x))]


def f2(ind,A,x):
	return [np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind])),np.fft.fft(A.dot(x[:,ind]))]

if __name__ == '__main__':
    pool = mp.Pool()
    print(mp.cpu_count())
    print(pool._processes)
    B  = np.random.rand(1000,1000)
    x0 = np.random.rand(1000,64) 
    indices = range(64)
    
    partial_f = partial(f1,A=B)

    a = time.time()
    xs = [] 
    for i in range(64):
        xs.append(x0[:,i]) 

    pool.map(partial_f,xs)
    b = time.time()
    print(b-a)

    a = time.time()
    for i in indices:
	f2(i,B,x0)
    b = time.time()
    print(b-a)

    pool.close()
    pool.join()


