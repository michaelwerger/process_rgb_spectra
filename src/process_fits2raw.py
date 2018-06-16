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
    if os.path.exists(rawfile):
        with open(rawfile, 'rb') as f:
            data = pickle.load(f)
            raw_r = data['raw_r']
            raw_g = data['raw_g']
            raw_b = data['raw_b']
    else:
            
        raw_r = np.zeros(maxcol)
        raw_g = np.zeros(maxcol)
        raw_b = np.zeros(maxcol)


        for data, fname in ic.data( return_fname=True):
            minrow = int(config.get('MAIN','ccdrows_min'))
            maxrow = int(config.get('MAIN','ccdrows_max'))

            rgb = demosaicing_CFA_Bayer_bilinear(data)
            
            trace = rgb[minrow:maxrow,mincol:maxcol,1].sum(axis=1)
            max_i = max(trace)
            ix_line = np.where(trace > (max_i*0.1))
            minrow = np.min(ix_line)
            maxrow = np.max(ix_line)

            print (fname, minrow, maxrow, maxrow-minrow, max_i)
            raw_r = raw_r + rgb[minrow:maxrow,mincol:maxcol,0].sum(axis=0)
            raw_g = raw_g + rgb[minrow:maxrow,mincol:maxcol,1].sum(axis=0)
            raw_b = raw_b + rgb[minrow:maxrow,mincol:maxcol,2].sum(axis=0)
        data = {
            'raw_r': raw_r,
            'raw_g': raw_g,
            'raw_b': raw_b,
        }
        with open(rawfile,'wb') as f:
            pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

    if args.dispersion == False and config.has_option(section,'coeff'):
        print ("Using existing dispersion relation")
        coeff = eval(config.get(section,'coeff'))
        einsen = np.linspace(1, len(raw_r), num=len(raw_r))
        wave = einsen * einsen * coeff[0] + einsen * coeff[1] + coeff[2]
        wave = wave[0:len(einsen)]

        fig = plt.figure(figsize=(14,6))
        plt.plot(wave,raw_r[::-1],color='r')
        plt.plot(wave,raw_g[::-1],color='g')
        plt.plot(wave,raw_b[::-1],color='b')

        for wave0 in eval(config.get(section,'lines')):
            plt.plot([wave0,wave0],[0.1e7,1.0e7])
            plt.text(wave0,1.0e7,str(wave0),rotation=90,rotation_mode='anchor', fontsize=10)
            #plt.Text(wave0,1.0e7,text=str(wave0))
        plt.xlim(3500,8000)
        plt.show()

        fig = plt.figure(figsize=(14,6))
        plt.plot(wave,raw_r[::-1]+2*raw_g[::-1]+raw_b[::-1])
        plt.show()
    else:
        print ("Mark spectrum lines for dispersion solution")
        print ("Hint: typical lines:")
        print ("Balmer Series: H-alpha 6563 H-beta 4861 H-gamma 4340 H-delta 4102 H-epsilon 3970 H-zeta 3889 H-eta 3835 ")
        print ("Telluric lines: 6863, 7594")
        fig = plt.figure(figsize=(14,6))
        plt.axis([1000, 3500, 0, 2.0e7])
        n = len(raw_r[::-1])
        pixels = np.linspace(1,n-1,num=n-1)

        pixels = np.arange(1,n+1,1)
        line, = plt.plot(pixels,raw_r[::-1],color='r')
        dc = DispersionDataCursor(plt.gca())
        fig.canvas.mpl_connect('pick_event', dc)
        line.set_picker(5) # Tolerance in points
        plt.show()

        print ("detected:")
        print (dc.getPositions())
        print (dc.getWavelengths())
        sortedPositions = np.sort(dc.getPositions())
        sortedWavelengths = np.sort(dc.getWavelengths())
        print (sortedPositions)
        print (sortedWavelengths)
        
        query = input('accept (Y/N)? ')
        if query == 'Y':
            config.set(section,'positions',value=str(sortedPositions))
            config.set(section,'wavelengths',value=str(sortedWavelengths))
            fig = plt.figure(figsize=(14,6))
            plt.plot(sortedPositions, sortedWavelengths)
            plt.show() 
            with open(configfile,'w') as cf:
                config.write(cf)


if __name__ == "__main__":
    main()