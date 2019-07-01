#NAME:
#   SPIDERS
#PURPOSE:
#   Model the spiders that support a Cassegrain secondary mirror for
#   purposes of WFE/PSF/MTF estimation. The resulting mask is to be
#   ADDED to a WFE map in order to block the wavefront. Note that
#   the mirror procedure is capable of setting up a circular aperture
#   with a central obstruction. Thus, MIRROR and SPIDERS together form
#   a comprehensive way of modeling the aperture of a Cassegrain telescope.
#CALLING SEQUENCE:
#   mask = spiders(N,widths,angles)
#EXAMPLE:
#   wfe += spiders(N,widths,angles) #the mask is ADDED to the WFE.
#INPUT PARAMETERS:
#   N = size of square WFE array (N x N)
#   widths = Array of spider widths in pixels.
#   angles = Array of spider angles in degrees. Must have same number
#       of elements as the widths array. Angles are measured CCW
#       from a horizontal axis drawn from mask array center to right edge.
#OUTPUT PARAMETERS:
#   mask = An NxN complex array whose real part is 1.0 everywhere. The
#       imaginary part is 1e6 or greater on the spiders and 0.0 elsewhere. When
#       mask is multiplied by a WFE map, it has the effect of zeroing the
#       electric field on the spiders.
#MODIFICATION HISTORY:
#   2008-Nov-14  C. Kankelborg
#   2019 Jun 21  JTE translated into Python 3

import numpy as np
from scipy.ndimage import rotate

def spiders(N,
            widths,
            angles):
    
    Nspiders = widths.size
    
    if angles.size != Nspiders:
        print('Number of widths =/= number of angles.')
    
    #Initialize the mask.
    mask = np.zeros([N, N], dtype = np.complex)
    
    #Create each spider and add it into the mask.
    for i in range(Nspiders):
        
        spider = np.zeros([N, N], dtype = np.complex)
        
        baseline = np.round(np.divide(N-widths[i] - 0.0,
                                      2))
        
        spider[N/2:N, baseline:round(baseline+widths[i])] = 1e6j
        
        mask += rotate(spider, -angles[i])
    
    return mask