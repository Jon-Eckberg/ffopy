#NAME:
#  mtf1d.pro
#PURPOSE:
#  Calculate modulation transfer function based on 1D line spread data. Take 
#  into account the profile of the source.
#CALLING SEQUENCE:
#  mtf = mtf1d(source, image)
#OUTPUTS:
#  mtf = modulation transfer function (MTF). Results are independent of the
#     normalization of source and image (even if they are different), because
#     the MTF is normalized to unit DC response. Of course, that's a limitation
#     of MTF as a measure of system performance.
#INPUTS:
#  source = n-element 1D array of intensities representing an ideal image of the
#     stimulus. In the simplest case, this would be an impulse function (one hot
#     pixel, the rest zeroes).
#  image = n-element 1D array of intensities representing how the optical system
#     images the source. If source is an impulse, then image is the
#     linespread function.
#OPTIONAL KEYWORD INPUTS:
#  dx = Pixel size. Default = 1 (i.e., distance units are pixels by default).
#  x_units = String describing the units for dx. Default='pixel' if dx
#     is not specified If dx is provided and x_units is not, then x_units
#     is ''. If x_units is defined, then it is used to infer k_units.
#  plotmtf = If set, then produce a plot of the MTF.
#  title = Passed straight through to IDL's plot procedure.
#  subtitle = Passed to IDL's plot procedure.
#  denoise = if set, assume the image signal at the Nyquist frequency
#     represents the noise floor, and subtract that level from the spectrum.
#     Denoise is set to the number of frequencies (including the Nyquist) to
#     be considered noise, and the mean amplitude of those is what's subtracted.
#OPTIONAL KEYWORD OUTPUTS:
#  ks = wavenumber axis for mtf, in cycles per distance unit. The distance units
#     are the same as those used for dx (pixels if dx is not supplied).
#  k_units = string describing the units of kx, ky, and ks. Constructed from
#     x_units, if it exists. If the units cannot be inferred, then k_units
#     is set to an empty string.
#  noise_floor = dimensionless noise floor (relative to MTF) as a function
#     of ks.
#MODIFICATION HISTORY:
#  2012-Jul-19  C. Kankelborg
#  2012-Aug-01  C. Kankelborg added noise_floor keyword.
#  2019 Jun 20  JTE translated to Python 3

import numpy as np
import matplotlib.pyplot as plt

def mtf1d(source,
          image,
          dx = 1.0,
          x_units = 'pixel',
          plotmtf = False,
          title = 'MTF1D',
          subtitle = '',
          denoise = False):

    k_units = 'cycles/{}'.format(x_units)
    
    n = source.size()
    
    if image.size() != n: print('Sizes of source:{0} and image:{1} are not equal.'.format(n, image.size()))
    
    
    source_f = np.abs(np.fft.fft(source))
    image_f  = np.abs(np.fft.fft(image))
    
    mtf = np.abs(np.divide(image_f,
                           source_f))
    
    #normalization
    
    mtf /= mtf[0]
    
    #Suppress noise that dominates at high frequencies
    if denoise:
        
       nyquist = np.floor(0.5 * n) #index of highest frequency
    
       hi_frequencies = nyquist - np.arange(denoise)
    
          #range of frequencies to average to find noise floor
       noise_floor_power = np.square(image_f[hi_frequencies]).mean()
    
          #mean *power* of noise floor
       image_f = np.sqrt((np.square(image_f) - noise_floor_power) > 0)
        
       #We're now in a position to calculate the noise floor amplitude (1 sigma), 
       #rescaled for comparison to MTF, as function of ks.
       noise_floor = np.divide(mtf[0] * np.sqrt(noise_floor_power),
                               source_f)
       
    
    #wavenumber axis in units of cycles per distance, where distance
    #is in the same units as dx.
       
    xshift = np.ceil(0.5 * n) - 1
    
    ks = np.divide(np.arange(n) - xshift,
                   n * dx) #Max is Nyquist frequency, 1/(2*dx).
    
    ks = np.shift(ks, -xshift)
    
    
    #Get rid of negative frequencies
    ss = np.where(ks >= 0)
    
    mtf = mtf[ss]
    
    ks = ks[ss]
    
    noise_floor = noise_floor[ss]
    
    if plotmtf:
        
        #xtitle = k_units
    
        #ytitle = 'MTF'
    
        plt.plot(ks, mtf, ylim=[0,1], label = 'MTF')
        
        if denoise: #plot noise floor as a dotted line
        
            plt.plot(ks, noise_floor, '.', label = 'noise floor')
    
        plt.legend()
    
    return mtf
    