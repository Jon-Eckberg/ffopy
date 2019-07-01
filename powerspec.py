#NAME:
#  powerspec
#PURPOSE:
#  Compute power spectrum of a 1D array by dividing into multiple overlapping segments
#  and averaging the power spectra. This improves the estimate at the expense of
#  spectral resolution. The normalization convention is the same as with 
#  periodogram_fft.pro. The method is a particular implementation of Welch's method
#  (IEEE Trans. Audio & Electroacoustics, 15:2:70-73, June, 1967). I have used the
#  Hanning cosinusoidal window.
#CALLING SEQUENCE:
#  result = powerspec(data [, Nfreq])
#INPUT PARAMETERS:
#  data = 1d array.
#  Nfreq = Number of frequencies (both positive and negative) to be used. Default=16, 
#     min=4.
#OPTIONAL KEYWORD INPUTS:
#  dx = sample interval. Default=1. If included, then k will be converted into units
#     of cycles (or radians, if that keyword is set) per unit distance (or time),
#     corresponding to the units of dx.
#  radians = if set, then let k be in radians (rather than cycles) per whatever.
#  retain_dc = if set, then retain the DC offset when calculating the power spectrum.
#     Ordinarily, the DC offset is subtracted before calculating the spectrum, because
#     the Hanning window tends to make the offset bleed into the nearest neighbors of
#     pixel [0,0] of the 2D power spectrum.
#OPTIONAL KEYWORD OUTPUTS:
#  kx = Nx-dimensional array of horizontal wavenumbers
#HISTORY:
#  2013-Aug-20 CCK
#  2013-Oct-03 CCK Added checking for infinities. A segment with any infinities is
#     excluded from the average power spectrum.


import numpy as np

def powerspec(data,
              Nfreq = 16,
              dx = dx,
              radians = False,
              retain_dc = False):
    
    N = data.size()
    
    if Nfreq < 4:
        raise ValueError('Nfreq is less than 4. Setting Nfreq = 4')
        Nfreq = 4 #Don't allow non-integer or smaller than 4.
    
    Nseg = np.ceil(np.divide( 2.0 * N, Nfreq )) - 1
       #How many segments  of length Nfreq fit in the interval, if we overlap them like
       #2 courses of bricks?
    if Nseg < 2:
        print('Data array too short to compute power spectrum using specified segment size.')
        
    xsep = np.divide(N - Nfreq, Nseg-1) * (1 - 1e-6) 
       #Except for the last factor, this should be the exact spacing required to fit.
       #The last factor is to ensure we can fit Nseg segments within the data array.
    spectra = np.full(Nseg, Nfreq)
    
    for i in range(Nseg):
        
       x0 = np.round(i*xsep)        #lower index of the segment
       x1 = x0 + Nfreq - 1  #upper index of the segment
       segment = data[x0:x1]
       spectra[i] = periodogram_fft(segment, retain_dc = retain_dc)
    
    #If there's a NaN or an Inf somewhere, we need to ignore it.
    ss = np.isfinite(spectra[:, 0]) #locations of valid (finite) spectra
    
    Nvalid = ss.size()
        
    power = np.divide(spectra[ss].sum(axis = 1), Nvalid)

    kx = k_arr(Nfreq,
               dx = dx,
               radians = radians)
    
    return power#, k_arr