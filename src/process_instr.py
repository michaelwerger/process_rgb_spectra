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
from InstrDataCursor import InstrDataCursor
import argparse
from scipy import interpolate

def main():


    configfile = 'processing.cfg'
    parser = argparse.ArgumentParser(description='Process raw RGB FITS files for specific star')
    parser.add_argument('--target', dest='target', default='ALPLYR', help='ALPLYR|GAMCYG')
    parser.add_argument('--configfile', dest='configfile', default=configfile, help='name for config file')
    
    
    args = parser.parse_args()

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
        print ("no input data")
        exit()
    
    instr_rgb_file = config.get(section,'instr_rgb_file')
    if os.path.exists(instr_rgb_file):
        with open(instr_rgb_file, 'rb') as f:
            instr_rgb = pickle.load(f)
            wave     = instr_rgb['wave']
            intens_r = instr_rgb['instr_r']
            intens_g = instr_rgb['instr_g']
            intens_b = instr_rgb['instr_b']

    if config.has_option(section,'coeff'):
        print ("Plotting raw data using existing dispersion relation")
        coeff = eval(config.get(section,'coeff'))
        einsen = np.linspace(1, len(raw_r), num=len(raw_r))
        wave = einsen * einsen * coeff[0] + einsen * coeff[1] + coeff[2]
        wave = wave[0:len(einsen)]

        fig = plt.figure(figsize=(14,6))
        line, = plt.plot(wave,raw_r[::-1],color='r')
        plt.plot(wave,intens_r,color='y')
        dc = InstrDataCursor(plt.gca())
        fig.canvas.mpl_connect('pick_event', dc)
        line.set_picker(5) # Tolerance in points
        plt.xlim(3500,8000)
        plt.ylim(0,max(raw_r[::1]*1.3))
        plt.show()
        print (dc.getPositions())
        print (dc.getIntensities())    
        if len(dc.getPositions()) > 1:
            tck = interpolate.splrep(dc.getPositions(),dc.getIntensities(),s=0)
            intens_r = interpolate.splev(wave,tck,der=0)
            fig = plt.figure(figsize=(14,6))
            plt.plot(wave,raw_r[::-1],color='r')
            plt.plot(wave,intens_r,'k-')
            plt.xlim(3500,8000)
            plt.ylim(0,max(raw_r[::1]*1.3))
            plt.show()
        
        fig = plt.figure(figsize=(14,6))
        line, = plt.plot(wave,raw_g[::-1],color='g')
        plt.plot(wave,intens_g,color='y')
        dc = InstrDataCursor(plt.gca())
        fig.canvas.mpl_connect('pick_event', dc)
        line.set_picker(5) # Tolerance in points
        plt.xlim(3500,8000)
        plt.ylim(0,max(raw_g[::1]*1.3))
        plt.show()
        print (dc.getPositions())
        print (dc.getIntensities())
        if len(dc.getPositions()) > 1:
            tck = interpolate.splrep(dc.getPositions(),dc.getIntensities(),s=0)
            intens_g = interpolate.splev(wave,tck,der=0)
            fig = plt.figure(figsize=(14,6))
            plt.plot(wave,raw_g[::-1],color='r')
            plt.plot(wave,intens_g,'k-')
            plt.xlim(3500,8000)
            plt.ylim(0,max(raw_g[::1]*1.3))
            plt.show()
        
        fig = plt.figure(figsize=(14,6))
        line, = plt.plot(wave,raw_b[::-1],color='b')
        plt.plot(wave,intens_b,color='y')
        dc = InstrDataCursor(plt.gca())
        fig.canvas.mpl_connect('pick_event', dc)
        line.set_picker(5) # Tolerance in points
        plt.xlim(3500,8000)
        plt.ylim(0,max(raw_b[::1]*1.3))
        plt.show()
        print (dc.getPositions())
        print (dc.getIntensities())    
        if len(dc.getPositions()) > 1:
            tck = interpolate.splrep(dc.getPositions(),dc.getIntensities(),s=0)
            intens_b = interpolate.splev(wave,tck,der=0)
            fig = plt.figure(figsize=(14,6))
            plt.plot(wave,raw_b[::-1],color='b')
            plt.plot(wave,intens_b,'k-')
            plt.xlim(3500,8000)
            plt.ylim(0,max(raw_b[::1]*1.3))
            plt.show()

        query = input('accept (Y/N)? ')
        if query == 'Y':
            
            instr_rgb = {
                'wave': wave,
                'instr_r': intens_r,
                'instr_g': intens_g,
                'instr_b': intens_b,
            }
            with open(instr_rgb_file,'wb') as f:
                pickle.dump(instr_rgb, f, pickle.HIGHEST_PROTOCOL)


    else:
        raise NotImplementedError('Only instr for intensity-wavelength plots')


if __name__ == "__main__":
    main()