# coding: utf-8

from __future__ import division, print_function
from ccdproc import ImageFileCollection
from glob import glob
import sys
from scipy import stats
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

import matplotlib.widgets as widgets

def inverse_poly1d(coeff,y):
    
    # assuming a polynomial given by coeff of the following order
    # y = cofff[2] + coeff[1] * x + coeff[2] * x^2
    # y = a0 + a1 x + a2 x^2
    # (y - a0) / a2 = a1/a2 x + x^2
    # x^2 + 2 a1/(2 a2) x + a1^2/(2 a2)^2 - a1^2/(2 a2)^2 = (y - a0) / a2
    # (x + a1/(2 a2))^2 = (y - a0) / a2 + a1^2/(2 a2)^2
    # x + a1/(2 a2) = +/- sqrt( (y - a0) / a2 + a1^2/(2 a2)^2 )
    # x_1 = - a1/(2 a2) + sqrt( (y - a0) / a2 + a1^2/(2 a2)^2 )
    # x_2 = - a1/(2 a2) - sqrt( (y - a0) / a2 + a1^2/(2 a2)^2 )

    const1 = coeff[1]/(2 * coeff[0])
    const2 = const1**2

    print ('coeff:' , coeff)
    if len(coeff) == 3:
        #a2 = coeff[0]
        #a1 = coeff[1]
        #a0 = coeff[2]

        x_1 = - const1 + np.sqrt( (y - coeff[2]) / coeff[0] + const2 )
        x_2 = - const1 - np.sqrt( (y - coeff[2]) / coeff[0] + const2 )
    else:
        raise NotImplementedError('only polynomials of 2nd order are supported')
    return [x_1, x_2]


positions = [1292.6267281105988, 1343.3179723502303, 1391.705069124424, 1465.4377880184329, 1599.0783410138247, 1884.7926267281105, 2792.6267281105993, 2953.9170506912442, 3329.4930875576038]
wavelengths = [3798, 3888, 3970, 4102, 4340, 4861, 6563, 6863, 7594]
positions = [1345, 1467, 1599, 1884, 2792, 2953, 3331]
wavelengths = [3889, 4102, 4340, 4861, 6563, 6863, 7594]
fig = plt.figure(figsize=(14,6))
plt.plot(positions, wavelengths,'o')


slope, intercept, r_value, p_value, std_err = stats.linregress(positions,wavelengths)

print (slope, intercept)
# y = slope * x + intercept

y0 = slope * positions[0] + intercept
yn = slope * positions[-1] + intercept
plt.plot([positions[0],positions[-1]], [y0,yn])

z = np.polyfit(positions, wavelengths, 2)
p = np.poly1d(z)

plt.plot(positions,p(positions),'r')
plt.show()




print ( positions )
print ( p(positions) )
print ( np.polyval(p, positions))
x_1, x_2 = inverse_poly1d(z, p(positions))
print ( x_1 )
#print ( x_2)