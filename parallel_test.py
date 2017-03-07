# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 17:26:18 2017

@author: dbanco02
"""
import numpy as np
import multiprocessing as mp

def f(x):
    print(x)
    return x*x

if __name__ == '__main__':
    pool = mp.Pool()
    print(mp.cpu_count())
    print(pool.processes)
    x = range(0,20) 
    pool.map(f,x)