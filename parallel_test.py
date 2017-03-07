# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 17:26:18 2017

@author: dbanco02
"""
import numpy as np
from multiprocessing import Pool

def f(x):
    return x*x

if __name__ == '__main__':
    pool = Pool(processes=6)
    x = [1,2,3,4,5,6,7,8]
    print(pool.map(f,x))