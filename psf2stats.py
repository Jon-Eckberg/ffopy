#NAME:
#  psf2stats
#PURPOSE:
#  Convert a point spread function to RMS spot radius (and optionally some
#  other statistics).
#CALLING SEQUENCE:
#  rmsrad = psf2stats(psf [, dx=dx])
#INPUT PARAMETERS:
#  psf = point spread function (floating point 2D array, Nx x Ny). This
#     need not be properly normalized.
#OPTIONAL KEYWORD INPUTS:
#  dx = PSF pixel size, in whatever length units are desired. Default=1.
#     Square pixels are assumed.
#  debug = If set, then break into immediate mode just before return.
#OPTIONAL KEYWORD OUTPUTS:
#  I0 = zeroth moment (total fluence just total(psf)).
#  xc, yc = centroid coordinates (origin at image center).
##  compactness1[DEPRECATED] = Experimental measure of PSF compactness. Bigger is more
##     compact. My intent here is to come up with a measure that is more
##     robust than RMS radius for analyzing focus data.
#  equivalent_area = PSF area in pixel^2. This is much more robust than
#     RMS radius for tasks like analyzing focus data. Conceptually analogous to equivalent 
#     width of a spectral line: https://en.wikipedia.org/wiki/Equivalent_width
#MODIFICATION HISTORY:
#  2009-Jun-23  C. Kankelborg
#  2016-Dec-22  CCK implemented dx default. Added compactness1.
#  2016-Dec-30  CCK replaced compactness1 with EQUIVALENT_AREA.
#  2019 Jun 20  JTE translated into Python 3

import numpy as np

def psf2stats(psf,
              dx = 1.0,
              I0 = I0,
              xc = xc,
              yc = yc,
              debug = False):
 
    Nx, Ny = psf.shape()
    
    #Coordinates (origin at image center):
    x = np.divide(np.arange(Nx) - (Nx-1),
                  2.0) * dx # replicate(1.0, Ny)
    y = np.ones(Nx) # (findgen(Ny) - (Ny-1)/2.0)*dx
    
    #Zeroth moment (total fluence):
    I0 = psf.sum()
    
    #First moments (centroid):
    xc = np.divide(x * psf, I0).sum()
    yc = np.divide(y * psf, I0).sum()
    
    #Second moments:
    r2 = np.square([x - xc, y - yc]).sum() #square of radius from spot center
    
    rmsrad = np.sqrt(np.divide(r2*psf,
                               I0).sum() )
    
    #Robust measure of spot size:
    smoothed = ck_smooth(psf, 1) 
       #boxcar smoothing alleviates inconsistencies due to sampling.
    #compactness1 = -alog(max(smoothed)/total(smoothed))
    equivalent_area = np.divide(smoothed.sum(),
                                smoothed.max())
    
    if debug: print('DEBUG option engaged: going to immediate mode.')
    
    return rmsrad#, equivalent_area
 