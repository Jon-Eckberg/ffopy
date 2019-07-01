#NAME:
#  MTF2D
#PURPOSE:
#  This is a general-purpose program to calculate the 2D modulation transfer
#  function (MTF) from an image of the point spread function (PSF). This is
#  much more informative than a 1D MTF, which only applies to a single
#  (often unspecified) orientation of the input plane-wave.
#CALLING SEQUENCE:
#  mtf = mtf2d(psf [, kx=kx, ky=ky]
#OUTPUT:
#  MTF = 2D image of modulation transfer function, as a function of kx and
#     ky. MTF is shifted so that the origin is in the middle. The optional
#     outputs kx and ky specify the precise alignment of the MTF array.
#INPUT:
#  PSF = 2D image of the point spread function. Need not be normalized
#     nor centered.
#OPTIONAL KEYWORD INPUTS:
#  dx = Pixel size. Default = 1 (i.e., distance units are pixels by default).
#  x_units = String describing the units for dx. Default='pixel' if dx
#     is not specified If dx is provided and x_units is not, then x_units
#     is ''. If x_units is defined, then it is used to infer k_units.
#  PLOT = If set, then produce a contour plot of the MTF.
#  TITLE = Passed straight through to IDL's CONTOUR procedure.
#  SUBTITLE = Passed to IDL's CONTOUR procedure.
#OPTIONAL KEYWORD OUTPUTS:
#  kx,ky = Coordinates for the MTF image. Units are cycles per unit
#     distance, where the distance units match the input dx.
#  k = magnitude of the spatial frequency, provided for convenience.
#  k_units = string describing the units of kx, ky, and k. Constructed from
#     x_units, if it exists. If the units cannot be inferred, then k_units
#     is set to an empty string.
#MODIFICATION HISTORY:
#  2012-Jun-26  C. Kankelborg
#  2012-Jul-10 CCK Fixed bug where odd-sized arrays could
#     have incorrect k axis. e.g.: xshift = ceil(nxp/2)-1
#     should be xshift = ceil(nxp/2.0)-1

import numpy as np
import matplotlib.pyplot as plt

def mtf2d(psf,
          dx = 1.0,
          plot = False,
          kx = kx,
          ky = ky,
          k = k,
          title = title,
          subtitle = subtitle,
          x_units = 'pixel'):

    k_units = 'cycles/{}'.format(x_units)
    
    #Get dimensions of the PSF array.
    nxp, nyp = psf.shape()
    
    mtf = np.abs(np.fft.fft2(psf))
    
    mtf /= mtf[0,0]
    
    xshift = np.ceil(nxp/2.0)-1
    
    yshift = np.ceil(nyp/2.0)-1
    
    mtf = np.shift(mtf, xshift, yshift) #2D MTF
    
    #wavenumber axes in units of cycles per distance, where distance
    #is in the same units as dx.
    kx = (np.arange(nxp) - xshift)/(nxp*dx) #Max is Nyquist frequency, 1/(2*dx).
    ky = (np.arange(nyp) - yshift)/(nyp*dx) #Max is Nyquist frequency, 1/(2*dx).
    k = np.sqrt(np.square([kx, ky]).sum())
    
    if plot:
       xtitle = 'kx {}'.format(k_units)
       ytitle = 'ky {}'.format(k_units)
       plt.contour(kx, ky, mtf, levels = [0.05, 0.1, 0.2, 0.4, 0.8])
       #xtitle=xtitle, ytitle=ytitle, title=title, subtitle=subtitle)
        
    return mtf
