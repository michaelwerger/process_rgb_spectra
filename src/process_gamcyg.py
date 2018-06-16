# coding: utf-8

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
from DataCursor import DataCursor


plt.style.use(astropy_mpl_style)

config = RawConfigParser()
config.read('processing.cfg')

section = 'GAMCYG'
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

coeff = eval(config.get(section,'coeff'))
einsen = np.linspace(1, len(raw_r), num=len(raw_r))
wave = einsen * einsen * coeff[0] + einsen * coeff[1] + coeff[2]
wave = wave[0:len(einsen)]
print( len(einsen))
print ( len(wave))

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

# SnapToCursor

#fig, ax = plt.subplots()

###cursor = Cursor(ax)
#cursor = SnaptoCursor(ax, wave,raw_r[::-1])
#cid =  plt.connect('motion_notify_event', cursor.mouse_move)

#ax.plot(wave,raw_r[::-1],)
#plt.axis([3500, 8000, 0, 2.0e7])
#plt.show()


# DataCursor

interactive = eval(config.get(section,'interactive'))

if interactive:
    print ("Mark spectrum lines")
    fig = plt.figure(figsize=(14,6))
    plt.axis([1000, 3500, 0, 2.0e7])
    n = len(raw_r[::-1])
    pixels = np.linspace(1,n-1,num=n-1)

    pixels = np.arange(1,n+1,1)
    line, = plt.plot(pixels,raw_r[::-1],color='r')
    dc = DataCursor(plt.gca())
    fig.canvas.mpl_connect('pick_event', dc)
    line.set_picker(5) # Tolerance in points
    plt.show()

    print ("detected:")
    print (dc.getPositions())
    print (dc.getWavelengths())

    fig = plt.figure(figsize=(14,6))
    plt.plot(dc.getPositions(), dc.getWavelengths())
    plt.show()
else:
    positions= eval(config.get(section,'positions'))
    wavelengths = eval(config.get(section,'wavelengths'))

    fig = plt.figure(figsize=(14,6))
    plt.plot(positions, wavelengths)
    plt.show()



