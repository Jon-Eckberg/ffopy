#NAME:
#  powerspec2d
#PURPOSE:
#  Compute power spectrum of a 2D array by dividing into many overlapping sub-images
#  and averaging the power spectra. This improves the estimate at the expense of
#  spectral resolution. The normalization convention is the same as with 
#  periodogram_fft.pro. This is a 2D analog of powerspec.pro, which uses Welch's
#  method with a Hanning window. Any sub-image containing NaN is ignored by
#  virtue of the /NAN keyword to mean() in periodogram_fft.pro.
#CALLING SEQUENCE:
#  result = powerspec2d(data [, Nfreq])
#INPUT PARAMETERS:
#  data = 2d array.
#  Nfreq = Size of square sub-image to analyze (optional). Default=16, min=4.
#OPTIONAL KEYWORD INPUTS:
#  dx = sample interval. Default=1. If included, then k will be converted into units
#     of cycles (or radians, if that keyword is set) per unit distance (or time),
#     corresponding to the units of dx.
#  dy = sample interval along the vertical direction (i.e., along the second
#     array index). Default: dy=dx.
#  radians = if set, then let k be in radians (rather than cycles).
#  retain_dc = if set, then retain the DC offset when calculating the power spectrum.
#     Ordinarily, the DC offset is subtracted before calculating the spectrum, because
#     the Hanning window tends to make the offset bleed into the nearest neighbors of
#     pixel [0,0] of the 2D power spectrum.
#OPTIONAL KEYWORD OUTPUTS:
#  kx = Nx-dimensional array of horizontal wavenumbers
#  ky = Ny-dimensional array of vertical wavenumbers
#  big_kx, big_ky = Nx x Ny arrays of horizontal and vertical wavenumbers, respectively.
#     This is handy for constructing digital filters as a function of the vector
#     wavenumber, to be applied to the FFT of an image.
#  k_magnitude = 2d array of wavenumber magnitudes.
#HISTORY:
#  2013-Aug-20 CCK
#  2013-Aug-21 CCK fixed dc offset issue, and added retain_dc keyword.
#  2015-Feb-21 CCK added /NaN in calls to mean(). I have verified that
#     the normalization is not affected by this change.
#  2019 Jun 25 JTE translated into Python 3

import numpy as np

def powerspec2d(data,
                Nfreq = 16,
                retain_dc = False,
                dx = dx,
                dy = dy,
                kx = kx,
                ky = ky,
                big_kx = big_kx,
                big_ky = big_ky,
                radians = radians):
    
    Nx, Ny = data.shape()
  
    Nfreq = round(Nfreq[0]) > 4 #Don't allow non-integer or smaller than 4.
    
    NsegX = np.ceil(np.divide(2.0 * Nx, Nfreq)) - 1.
    
    NsegY = np.ceil(np.divide(2.0 * Ny, Nfreq)) - 1.
       #How many segments  of length Nfreq fit in the interval, if we overlap them like
       #2 courses of bricks?
    
    if (NsegX < 2 | NsegY < 2):
        print('Data array too short to compute power spectrum using specified segment size.')
    
    xsep = np.divide(Nx - Nfreq, NsegX - 1) * (1 - 1e-6) 
    ysep = np.divide(Ny - Nfreq, NsegY - 1) * (1 - 1e-6) 
       #The last factor is to ensure we can fit Nseg segments within the data array.
    
    spectra = np.full([NsegX, NsegY], Nfreq)
    
    for i in range(NsegX):
        
       x0 = np.round(i * xsep)   #lower x index of the segment
       
       x1 = x0 + Nfreq - 1  #upper x index of the segment
       
       for j in range(NsegY):
           
          y0 = np.round(j * ysep)   #lower y index of the segment
          
          y1 = y0 + Nfreq - 1  #upper y index of the segment
          
          segment = data[x0:x1, y0:y1]
          
          spectra[i, j] = periodogram_fft(segment,
                                          retain_dc = retain_dc)
    
    power = np.nanmean(np.nanmean(spectra, axis = 1), axis = 1)
    
    #pick off leading dimension twice!
    
    k_magnitude = k_arr2d(Nfreq,
                          Nfreq,
                          dx =dx,
                          dy =dy,
                          radians =radians,
                          kx =kx,
                          ky =ky,
                          big_kx =big_kx,
                          big_ky =big_ky)
      
    return power #, k_magnitude