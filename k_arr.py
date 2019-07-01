#NAME:
#  K_ARR
#PURPOSE:
#  Generate an array of wavenumber values (by default, in units of cycles per 
#  sample interval), arranged corresponding to the FFT of an N-element array.
#CALLING SEQUENCE:
#  k = k_arr(N [, dx=dx] [,/radians])
#OUTPUT:
#  k = N-element array of positive and negative wavenumbers, where the maximum
#     (Nyquist) frequency has a value of pi radians per sample interval. The
#     dx and cycles keywords can be used to modify the units of k.
#INPUTS:
#  N = Number of elements. The program has been verified for both odd and
#     even N.
#OPTIONAL KEYWORD INPUTS:
#  dx = sample interval. Default=1. If included, then k will be converted into units
#     of cycles (or radians, if that keyword is set) per unit distance (or time),
#     corresponding to the units of dx.
#  radians = if set, then let k be in radians (rather than cycles).
#HISTORY:
#  2013-Jun-25 C. Kankelborg
#  2019 Jun 21  JTE Translated into Python 3

import numpy as np

def k_arr(N,
          dx = 1.0,
          radians = False):
    
    dk = np.power(N * dx, -1.0)  #frequency interval equals the fundamental frequency
       #of one cycle per N samples, or one cycle per total distance N*dx.
    
    k = dk * [ np.arange(np.floor_divide(N, 2.0) + 1), -np.arange(np.ceil(0.5 * N - 1))[::-1] + 1  ]
    #        [  positive frequencies,        negative frequencies              ]
    
    if radians:
        
        k *= 2.0 * np.pi #convert from cycles to radians.
    
    return k
