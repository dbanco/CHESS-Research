# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:18:44 2024

@author: dpqb1
"""

# process_spot.py

import sys
import spotfetch as sf

def main(k, s, t, params, outPath, fname1, fname2):
    sf.processSpot(k, s, t, params, outPath, fname1, fname2)

if __name__ == "__main__":
    k = int(sys.argv[1])
    s = int(sys.argv[2])
    t = int(sys.argv[3])
    params = sys.argv[4]  # This might be a string that you need to parse back into a dictionary
    outPath = sys.argv[5]
    fname1 = sys.argv[6]
    fname2 = sys.argv[7]
    
    main(k, s, t, params, outPath, fname1, fname2)
