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
from scipy import interpolate
import csv

def main():

    '''
    F_target = I_target / s = I_target * F_source / instr

    '''

    configfile = 'processing.cfg'
    parser = argparse.ArgumentParser(description='Process raw RGB FITS files for specific star')
    parser.add_argument('--target', dest='target', default='GAMCYG', help='GAMCYG')
    parser.add_argument('--source', dest='source', default='ALPLYR', help='ALPLYR')
    parser.add_argument('--configfile', dest='configfile', default=configfile, help='name for config file')
    
    args = parser.parse_args()

    plt.style.use(astropy_mpl_style)

    config = RawConfigParser()
    config.read(configfile)

    target_section = args.target
    source_section = args.source
    fitsdirs = glob(config.get('MAIN','datapath'))
    
    instr_rgb_file = config.get(source_section,'instr_rgb_file')
    target_rawfile = config.get(target_section,'raw_rgb_file')

    mincol = int(config.get('MAIN','ccdcols_min'))
    maxcol = int(config.get('MAIN','ccdcols_max'))
    minrow = int(config.get('MAIN','ccdrows_min'))
    maxrow = int(config.get('MAIN','ccdrows_max'))

    print ("loading data ...")

    filename = '/Users/Micha/data/ref/hst/alphalyr.tab'
    w = []
    f = []
    startRow = 16
    with open(filename, ) as csvfile:
        fluxreader = csv.reader(csvfile, delimiter='\t')
        for i in range(0,startRow):
            fluxreader.__next__()
        for row in fluxreader:
            w.append(float(row[0][:10]))
            f.append(float(row[0][11:]))
            #print ( w, f)
    
    wavestart = 3500
    waveend = 8000
    wave_np = np.arange(wavestart,waveend,1)
    print (len(wave_np))
    print ( wave_np[0:9])
    print ( wave_np[-10:])

    w_np = np.array(w)
    f_np = np.array(f)
    ix1s = np.where((w_np[w_np < waveend] > wavestart))
    tck = interpolate.splrep(w_np[ix1s],f_np[ix1s],s=0)
    flux_np = interpolate.splev(wave_np,tck,der=0)
    
    fig = plt.figure(figsize=(14,6))
    plt.plot(w,f,'k')
    plt.plot(wave_np,flux_np,'r')
    plt.xlim(wavestart,waveend)
    plt.ylim(0,6e9)
    plt.show()

    
    if os.path.exists(instr_rgb_file):
        with open(instr_rgb_file, 'rb') as f:
            instr_rgb = pickle.load(f)
        print ("instr data  loaded")
        wave     = instr_rgb['wave']
        intens_r = instr_rgb['instr_r']
        intens_g = instr_rgb['instr_g']
        intens_b = instr_rgb['instr_b']

        tck = interpolate.splrep(wave,intens_r,s=0)
        instr_r = interpolate.splev(wave_np,tck,der=0)
        tck = interpolate.splrep(wave,intens_g,s=0)
        instr_g = interpolate.splev(wave_np,tck,der=0)
        tck = interpolate.splrep(wave,intens_b,s=0)
        instr_b = interpolate.splev(wave_np,tck,der=0)

        fig = plt.figure(figsize=(14,6))
        plt.plot(wave,intens_r,'r')
        plt.plot(wave,intens_g,'g')
        plt.plot(wave,intens_b,'b')
        plt.xlim(wavestart,waveend)
        plt.ylim(0,6e7)
        plt.show()
    else:
        print ("no source instr data")
        exit()

    if os.path.exists(target_rawfile):
        with open(target_rawfile, 'rb') as f:
            data = pickle.load(f)
        print ("target data  loaded")
        raw_r = data['raw_r'][::-1]
        raw_g = data['raw_g'][::-1]
        raw_b = data['raw_b'][::-1]

        tck = interpolate.splrep(wave,raw_r,s=0)
        meas_r = interpolate.splev(wave_np,tck,der=0)
        tck = interpolate.splrep(wave,raw_g,s=0)
        meas_g = interpolate.splev(wave_np,tck,der=0)
        tck = interpolate.splrep(wave,raw_b,s=0)
        meas_b = interpolate.splev(wave_np,tck,der=0)
        fig = plt.figure(figsize=(14,6))
        plt.plot(wave_np,meas_r,'r')
        plt.plot(wave_np,meas_g,'g')
        plt.plot(wave_np,meas_b,'b')
        plt.xlim(wavestart,waveend)
        plt.ylim(0,1e7)
        plt.show()
    else:
        print ("no target raw data")
        exit()

    min_i = 0
    max_i = 4500
    flux_r = np.linspace(min_i, max_i, num=max_i) * 0.0
    flux_g = flux_r
    flux_b = flux_r
    
    # try:
    #     for i in range(min_i, max_i):
    #         offset=1129
    #         j = i + offset
    #         flux_r[i] = flux_np[i] * meas_r[i]/instr_r[i]
    #         flux_g[i] = flux_np[i] * meas_g[i]/instr_g[i]
    #         flux_b[i] = flux_np[i] * meas_b[i]/instr_b[i]
            
    #         #print (i, j, wave_np[i], flux_np[i], intens_r[j], raw_r[j], flux_r[i])
    # except IndexError:
    #     pass    

    flux_r = flux_np * meas_r / instr_r
    flux_g = flux_np * meas_g / instr_g
    flux_b = flux_np * meas_b / instr_b
    
    flux_r = np.multiply(flux_np, meas_r)
    flux_r = np.divide(flux_r, instr_r)
    flux_g = np.multiply(flux_np, meas_g)
    flux_g = np.divide(flux_g, instr_g)
    flux_b = np.multiply(flux_np, meas_b)
    flux_b = np.divide(flux_b, instr_b)
    
    fig = plt.figure(figsize=(14,6))
    #plt.plot(wave_np,flux_r,'r')
    
    plt.plot(wave_np,flux_g,'g')

    #plt.plot(wave_np,flux_b,'b')
    plt.xlim(wavestart,waveend)
    plt.ylim(0,2e9)
    plt.show()


if __name__ == "__main__":
    main()