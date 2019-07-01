#NAME:
#  wfe2psf
#PURPOSE:
#  Convert a wavefront error map to a point spread function.
#CALLING SEQUENCE#
#  psf = wfe2psf(wfe [,oversample=oversample])
#       ...see below for more keyword options!
#INPUT PARAMETERS:
#  wfe = wavefront error in waves, which may be complex as documented 
#  in psd2wfe.pro. Must be a 2D square (NxN) array. Note that large
#  imaginary values essentially result in zero Efield (darkness). This
#  is how apertures are marked off. And by more subtle application,
#  apodization or nonuniform illumination effects could also be
#  implemented.
#OUTPUT PARAMETERS:
#  psf = (scalar, 2D array) point spread function, same dimensions
#     as the input wfe, and normalized to 100% response at DC.
#OPTIONAL KEYWORD INPUTS:
#  oversample = factor by which to oversample the PSF. Oversampling is
#     highly recommended. Default=4. Only integer values have been tested,
#     but maybe a non-integer would work. Do you feel lucky?
#  fratio = The f-ratio of the beam (focal length / beam width). The fratio
#     should be calculated using the total width of the wfe map, not
#     the size of a smaller aperture embedded in that map! Don't use
#     an effective f-ratio or T-ratio here. The fratio is used only for
#     geometric calculation of dx and xarr (see below). Note that small
#     angle approximation is used ragardless of fratio. If the value
#     1.0/width is put into the fratio parameter, where width is the width
#     of the wfe map in the same units as lam, then the dx and xarr
#     outputs will be in radians.
#  lam = Wavelength of the light. See dx and xarr below.
#OPTIONAL KEYWORD OUTPUTS:
#   dx = scale of pixels of the output image. Requires fratio and 
#       wavelength keywords. The units are the same as the wavelength units.
#       Alternatively, if fratio = 1.0/width is used, then dx
#       will be in radians (see fratio).
#   xarr = array of axis values for the side of the PSF image. Requires
#       fratio and wavelength keywords. Units same as wavelength.
#       Alternatively, if fratio = 1.0/width is used, then xarr will
#       be in radians (see fratio).
#MODIFICATION HISTORY:
#  2008-Nov-07  C. Kankelborg
#  2008-Nov-10  CCK added dx and xarr keywords.
#  2012-Nov-27  CCK clarified documentation a bit.
#  2015-Feb-14  CCK improved documentation. Clarified how to get PSF
#     coordinates in radians rather than distance units.
#  2019 Jun 20  JTE translated to Python 3

import numpy as np

def wfe2psf(wfe,
            oversample = 4,
            f_ratio = 0.0,
            lam = 0.0):
    
    N = wfe.size()
    
    bigmap = np.zeros([oversample*N,
                       oversample*N],
                      dtype = np.complex)
       #padded to the degree specified by the oversample parameter.
    bigmap[:N, :N] = np.exp(wfe * 2.0j * np.pi)
    
    bigmap = np.shift(bigmap, -N/2, -N/2)
    
    Efield = np.fft.fft2(bigmap)
    
    intensity = np.square(np.absolute(Efield))
    
    intensity_small = (np.shift(intensity, N/2, N/2) )[:N, :N]
    
    psf = np.divide(intensity_small,
                    intensity_small.sum())
        
    #Work out the scale of the PSF map.
    #Without oversampling, one pixel would subtend an angle equal to the
    #wavelength over the beam diameter.
    if (f_ratio > 0.0 & lam > 0.0):
        
        dx = np.divide(lam * f_ratio,
                       oversample)
        
        xarr = ( (N) - 0.5 * N) * dx
    
        return psf, dx, xarr 
    
    else:
        
        return psf 
