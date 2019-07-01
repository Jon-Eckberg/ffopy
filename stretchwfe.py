#NAME:
#  STRETCHWFE
#PURPOSE:
#  Stretch a WFE map while keeping the array size the same. 
#  The advantage of this approach is that the f/ratio of the WFE (ratio of 
#  focal length to WFE map width) remains the same, to avoid confusion. Use,
#  e.g., when beam geometry changes with a reflection off a diffraction grating.
#  If the stretching factor in either axis is greater than unity, then the 
#  aperture (beam size within WFE map) needs enough margin to move toward the 
#  edge without falling over. Use, for example, the OVERSIZE keyword to MIRROR 
#  or PSD2WFE.
#CALLING SEQUENCE:
#  wfe_new = stretchwfe(wfe, factor)
#  wfe_new = stretchwfe(wfe, [factor_x, factor_y])
#OUTPUT:
#  wfe_new -- Stretched WFE map.
#INPUTS:
#  wfe -- A WFE map.
#  factor -- 2D array of stretching factors in the horizontal and vertical 
#     directions. Greater than 1 stretches, less than 1 squeezes.
#MODIFICATION HISTORY:
#  2009-Oct-17  C. Kankelborg
#  2019 Jun 25  JTE translated into Python 3

from numpy import np
from scipy.ndimage import zoom

def stretchwfe(wfe,
               factor):
    
    Nx, Ny = wfe.shape
    
    #if min(factor) lt 1 then message,'FACTOR needs to be greater than 1. Exiting.'
    if factor.ndim != 2: print('FACTOR needs to be a 2d array. Exiting.')
    
    Nx_new = np.round(Nx * factor[0])
    Ny_new = np.round(Ny * factor[1])
    
    wfe_new = zoom(wfe, factor, order = 1) # Order 1 = linear
       # More pixels for same pattern
       #Linear interpolation is important here because the discontinuity at
       #the aperture edge would cause ringing if I used cubic interpolation.
    
    #if wfe shrinks, pad with opaque aperture
    if (Nx_new < Nx):
       foo = 1e9j * np.ones([np.divide(Nx - Nx_new + 1, 2), Ny_new])
       wfe_new = [foo, wfe_new, foo]
       Nx_new = wfe_new.shape[0]
    
    
    if Ny_new < Ny:
       foo = 1e9j * np.ones([Nx_new, np.divide(Ny - Ny_new + 1, 2)])
       wfe_new = [foo, wfe_new, foo]
       Ny_new = wfe_new.shape[1]
    
    #if wfe has expanded, cut it down. In fact, even if it shrunk, there may have
    #been a slight expansion afterward due to off-by-half in the above padding!
    x0 = 0.5 * (Nx_new - Nx)
    x1 = x0 + Nx - 1
    y0 = 0.5 * (Ny_new - Ny)
    y1 = y0 + Ny - 1
    
    return wfe_new[x0:x1, y0:y1]