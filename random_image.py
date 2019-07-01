#NAME:
#  RANDOM_IMAGE
#PURPOSE:
#  Generate a random image with a specified power law spectrum. The result is a floating
#  point, 2d array with all positive values. The result will be periodic unless padding
#  is used (see the pad keyword).
#CALLING SEQUENCE:
#  image = RANDOM_IMAGE(seed, Nx, Ny, alpha, meanval [, /complex])
#INPUT PARAMETERS:
#  seed = seed for the random number generator. See randomu.
#  Nx, Ny = size of the image (powers of 2 work best# code untested with odd values).
#  alpha = power law index for the image power spectrum (goes as k^(-alpha)).
#  meanval = mean value of the image.
#OPTIONAL KEYWORD INPUTS:
#  pad = positive integer factor for padding the arrays. pad=2 or higher prevents the
#     output image from being periodic in structure. If a non-integer value, or a value
#     less than 1 is input, then pad is set to the nearest integer greater than or equal
#     to 1.
#HISTORY:
#  2013-Jun-25 CCK 
#  2013-Jul-11 CCK revised documentation (seed argument had been omitted from call seq).
#  2013-Aug-20 CCK rewritten to produce image with a more accurately determined
#     power spectrum.
# 2019 Jun 21  JTE Translated into Python 3

import numpy as np

def random_image(seed,
                 Nx,
                 Ny,
                 alpha,
                 meanval,
                 pad = 1):
    
    k = k_arr2d(Nx * pad,
                Ny * pad,
                big_kx = kx,
                big_ky = ky) #2d array of scalar wavenumbers
    
    powerspec = np.power(k, -alpha) #power spectrum of the image
    powerspec[0, 0] = 0.0
    
    image = np.random.rand(Nx*pad, Ny*pad) 
    
    image_f = np.fft.fft2(image)
    
    image_power = np.square(np.absolute(image_f))
    
    image_f *= np.sqrt(np.divide(powerspec,
                                 image_power))
    
    image = np.fft.ifft2(image_f)[1:-1, 1:-1].real #index range removes padding.
    
    image -= image.min() #make it positive
    image *= meanval / image.mean() #adjust the mean to whatever was specified
    
    return image
