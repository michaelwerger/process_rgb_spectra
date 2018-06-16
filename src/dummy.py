#-*- coding: utf-8 -*-

from __future__ import division, print_function
from ccdproc import ImageFileCollection
from glob import glob
import sys
import numpy as np
import colour
from colour.plotting import *
import colour_demosaicing
from colour_demosaicing import (
    demosaicing_CFA_Bayer_bilinear,
    demosaicing_CFA_Bayer_Malvar2004,
    demosaicing_CFA_Bayer_Menon2007,
    mosaicing_CFA_Bayer)
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
import os
import pickle
from configparser import RawConfigParser

import matplotlib.widgets as widgets

from SnapToCursor import SnaptoCursor
from DispersionDataCursor import DispersionDataCursor
import argparse

def main():


    configfile = 'processing.cfg'
    parser = argparse.ArgumentParser(description='Process raw RGB FITS files for specific star')
    parser.add_argument('--target', dest='target', default='ALPLYR', help='ALPLYR|GAMCYG')
    parser.add_argument('--configfile', dest='configfile', default=configfile, help='name for config file')
    parser.add_argument('--dispersion', dest='dispersion', action='store_true', default=False, help='find dispersion relation by marking lines')
    
    args = parser.parse_args()

    print (args.dispersion)
    plt.style.use(astropy_mpl_style)

    config = RawConfigParser()
    config.read(configfile)

    section = args.target
    fitsdirs = glob(config.get(section,'datapath'))
    master = config.get(section,'master')
    numbers = eval(config.get(section,'numbers'))
    rawfile = config.get(section,'raw_rgb_file')
    mincol = int(config.get('MAIN','ccdcols_min'))
    maxcol = int(config.get('MAIN','ccdcols_max'))
    minrow = int(config.get('MAIN','ccdrows_min'))
    maxrow = int(config.get('MAIN','ccdrows_max'))

    filenames = []
    for number in numbers:
        filenames.append("%s%03d%s" % (master,number,'.fit'))

    ic = ImageFileCollection(fitsdirs[0], keywords='*', filenames=filenames)

    rawfile = config.get(section,'raw_rgb_file')
    
            
    raw_r = np.zeros(maxcol)
    raw_g = np.zeros(maxcol)
    raw_b = np.zeros(maxcol)


    for data, fname in ic.data( return_fname=True):
        minrow = int(config.get('MAIN','ccdrows_min'))
        maxrow = int(config.get('MAIN','ccdrows_max'))

        rgb = demosaicing_CFA_Bayer_bilinear(data)
        
        trace = rgb[minrow:maxrow,mincol:maxcol,1].sum(axis=1)
        max_i = max(trace)
        
        fig = plt.figure(figsize=(14,6))
        #plt.plot(wave_np,flux_r,'r')
        
        plt.plot(trace/max_i)

        
        plt.show()
    
if __name__ == "__main__":
    main()