#+
#NAME:
#  ASTIGMATISM
#PURPOSE:
#  Modify a complex wavefront error (WFE) image by adding a specific amount
#  of astigmatism.
#CALLING SEQUENCE:
#  wfe_new = astigmatism(wfe_orig, rms_astigmatism, clocking)
#INPUT PARAMETERS:
#  wfe_orig -- 2D complex image of initial WFE, in waves.
#  rms_astigmatism -- value of RMS WFE to be applied in the form
#     of astigmatism. If you want PV or something else, hey, scale it yourself.
#  clocking -- orientation in radians of the astigmatism pattern. Zero
#     corresponds to maximum values on the x-axis. Positive clocking is
#     counterclockwise according to the usual way of viewing images in IDL.
#OUTPUT PARAMETERS:
#  wfe_new -- wavefront error with the specified astigmatism added.
#OPTIONAL OUTPUT KEYWORDS:
#  wfe_change -- change in WFE (wfe_new - wfe_orig).
#  rmswfe -- RMS wavefront error of the result.
#ALGORITHM:
#  We calculate a wavefront error given by:
#     WFE_change =  r^2*(cos(theta-clocking)^2 - sin(theta-clocking)^2),
#  where r, theta are polar coordinates centered in the array.
#  This array is then scaled to the desired rms_astigmatism.
#  The magnitude of this scaling may vary based on the detailed aperture
#  map, which is specified by the imaginary pat of wfe_orig.
#BUGS:
#  None left, I hope.
#MODIFICATION HISTORY:
#  2009-Dec-10  C. Kankelborg
#  2010-Jan-11  CCK corrected documentation error: rms_astigmatism is
#     a WFE, not figure error!
#  2019 Jun 25   JTE translated into Python 3

import numpy as np

def astigmatism(wfe_orig,
                rms_astigmatism,
                clocking,
                rms_wfe = rms_wfe):

    #If clocking is a 1-element array, it will screw everything up. That's
    #because, when IDL multiplies a 1-element array by an n-element array,
    #the result is perversely a 1-element array (it's not treated as
    #scalar multiplication!). If, for example, clocking is the output of 
    #randomu, then this sort of thing will happen. It's stupid, but the
    #following line fixes the problem.
    clocking = clocking[0]
    
    #Create normalized aperture coordinates
    Nx, Ny = wfe_orig.size()
    
    xlin = np.divide(np.arange(Nx),
                     Nx - 1.0) - 0.5
    
    ylin = np.divide(np.arange(Ny) - 0.5 * (Ny - 1),
                     Nx - 1.0)
    
    x = xlin # replicate(1.0, Ny)
    
    y = np.ones(Nx) # ylin
    r = np.sqrt(np.square([x, y]).sum())
    
    theta = np.atan(y, x)
    
    #Calculate WFE for astigmatism, scaling as promised
    wfe_change = np.square(r) * np.square([np.cos(theta-clocking), -np.sin(theta-clocking)]).sum()
    
    powers = np.exp(-2.0 * wfe_orig.imag)
      
    powers /= powers.sum()
    
       #use aperture power map to specify weights for the RMS.
    
    rms_unscaled = np.sqrt(np.square(powers * wfe_change).sum())
    wfe_change *= np.divide(rms_astigmatism,
                            rms_unscaled)
    
    #Add astigmatism into the WFE
    wfe = wfe_orig + wfe_change
    
    rms_wfe = np.sqrt(np.square(powers * wfe.real).sum())
    
    return wfe#, rms_wfe