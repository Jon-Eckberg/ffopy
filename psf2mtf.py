#NAME:
#  psf2mtf
#PURPOSE:
#  Convert a point spread function to a modulation transfer function
#  (and optionally a line spread function).
#CALLING SEQUENCE:
#  mtf = psf2mtf(psf [, theta=theta] [, lsf = lsf] [, dx=dx], [, karr=karr])
#INPUT PARAMETERS:
#  psf = point spread function (floating point 2D array, Nx x Ny). This
#     need not be properly normalized.
#OUTPUT PARAMETERS:
#  mtf = modulation transfer function, which is the amplitude response
#     as a function of scalar wavenumber. This is a 1D array with Nx/2
#     elements. Normalized to unit DC response.
#OPTIONAL KEYWORD INPUTS:
#  theta = CCW orientation of sinusoidal test pattern in degrees. 
#     Default, theta=0, corresponds to vertical bars (horizontal 
#     resolution).
#  dx = PSF pixel size, in whatever length units are desired for karr.
#OPTIONAL KEYWORD OUTPUTS:
#  lsf = linespread function (which is just PSF summed along the orientation
#     specified by theta).
#  karr = frequency axis for the MTF
#     in periods per pixel (or per unit length if dx is supplied).
#MODIFICATION HISTORY:
#  2008-Nov-07  C. Kankelborg
#  2008-Nov-10  CCK added dx and karr keywords.
#  2015-Mar-30  CCK switched to using k_arr, standard method to
#     work out the wavenumbers. The karr output is now produced
#     even if dx is not given.
#  2019 Jun 21   JTE translated into Python 3


from scipy.ndimage import rotate
import numpy as np

def psf2mtf(psf,
            theta = 0,
            lsf = lsf,
            dx = dx,
            karr = karr):
    
    Nx = psf.shape[0]
    
    psf_rot = rotate(psf, theta)
    
    lsf = psf_rot.sum(axis = 2) #linespread function
    
    mtf = np.abs(np.fft.fft(lsf))[:Nx/2] #leave out the redundant negative frequencies
    
    karr = (k_arr(Nx,
                  dx = dx))[:Nx/2]
    
    return np.divide(mtf, mtf[0]) #properly normalized, unity DC response one hopes!