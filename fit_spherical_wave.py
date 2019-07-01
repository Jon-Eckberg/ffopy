#NAME:
#   fit_spherical_wave
#PURPOSE:
#   Find the best-fit spherical wave to a WFE map so that it can be subtracted 
#   out. This is equivalent to removing defocus, tip and tilt. Actually, I have
#   settled for a paraboloid approximation of a sphere. The best fit is defined
#   by minimizing the sum of the absolute deviations (robust fitting), unless
#   the program is invoked with the /least_squares or /linear option.
#CALLING SEQUENCE:
#   swave = fit_spherical_wave(wfe)
#INPUT PARAMETERS:
#   WFE = 2D wavefront error map. The units can be phase, waves, nanometers,
#   or whatever. Compatible with complex WFE as introduced in psd2wfe.
#OUTPUTS:
#   SWAVE = Best fit spherical surface, in the same units as WFE. The WFE
#       pixels are assumed to be square.
#OPTIONAL INPUT KEYWORDS:
#  linear --- If set, then use orthogonal polynomials. The catch is that,
#     if the aperture is too 'interesting', these may turn out not to be
#     orthogonal after all!
#  least_squares --- If set, then fit by least squares. This overrides the
#     default minimum absolute deviation approach. It would be much faster
#     to use linear least squares, but that was not my approach.
#     If you want fast, use /linear!
#MODIFICATION HISTORY:
#   2008-Nov-03  C. Kankelborg
#   2008-Nov-10  CCK adapted to complex wfe (see psd2wfe.pro).
#   2009-May-28  CCK tweaked precision and max iterations for 
#        use with moses_ideal.pro.
#   2011-Nov-07  CCK added /least_squares and /linear options.
#   2019 Jun 21   JTE translated to Python 3

import numpy as np

def paraboloid(K, x0, xp, y0, yp, z1):

    return 0.5 * K * np.square([xp-x0,
                                yp-y0]).sum() + z1


def badness(wfe,
            aperture,
            K,
            x0,
            xp,
            y0,
            yp,
            z1,
            least_squares = False):
    
    sphere = paraboloid(K, x0, xp, y0, yp, z1)
    
    if least_squares:
        
        result = np.square(sphere[aperture] - wfe[aperture]).sum()
    
    else:

        result = np.abs(sphere[aperture] - wfe[aperture]).sum() + 1.0 
          #The 1.0 aids in convergence!
    
    return result




def fit_spherical_wave(wfe,
                       wfe_aperture,
                       linear = True,
                       least_squares = False):

    aperture = np.where(wfe.imag < 1.) #adaptation to complex WFE
    
    Nx, Ny = wfe.shape
    
    #Create rectangular arrays for x and y coordinates.
    xp = np.outer(np.arange(Nx), np.ones(Ny))
    yp = np.outer(np.ones(Nx), np.arange(Ny))
    
    if linear:
       #Note that double precision is used for all the means over many
       #elements, to minimize roundoff error. The result, however, is
       #cast back to float before returning.
       
       result = wfe[aperture].mean() * np.ones((Nx, Ny), dtype = np.float) #piston term
       #Construct basis elements:
       b1 = xp.copy() #tip
       b1 -= b1[aperture].mean() #remove bias
       b1 /= np.sqrt(np.square(b1[aperture]).sum()) #normalize
       b2 = yp.copy() #tilt
       b2 -= b2[aperture].mean() #remove bias
       b2 /= np.sqrt(np.square(b2[aperture]).sum()) #normalize
       b3 = np.square(b1) + np.square(b2) #defocus
          #important to use b1 and b2 rather than xp and yp, so b3 not decentered!
          
       b3 -= b3[aperture].mean() #remove bias
       b3 /= np.sqrt(np.square(b3[aperture]).sum()) #normalize
       
       #Remove component along each basis element in turn:
       result += b1 * wfe_aperture * b1[aperture].sum()
       result += b2 * wfe_aperture * b2[aperture].sum()
       result += b3 * wfe_aperture * b3[aperture].sum()
       
       return result
    
    
    initial_prams = [0,0, 0.5 * Nx, 0.5 * Ny]
    
    scale_prams = [1.0, 1.0/((Nx/2)^2.0 + (Ny/2)^2.0), Nx, Ny]
    
    initial_prams += scale_prams * 0.1 * np.random.random(4)
       #slight perturbation -- CCK 20100407 -- trying to break pathologies
    
    #N-M simplex
    prams_fit = amoeba(1e-5,
                       FUNCTION_NAME = 'badness',
                       P0=initial_prams,
                       scale=scale_prams,
                       nmax=10000)
    
    return paraboloid(prams_fit)
    
    
    
