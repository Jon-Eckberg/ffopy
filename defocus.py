#+
#NAME:
#  DEFOCUS
#PURPOSE:
#  Modify a wavefront error (WFE) image to simulate an increase (or decrease)
#  in distance to the detector. Resulting PSFs obtained with WFE2PSF will
#  show the focus changes, and can be used to create a focus series.
#CALLING SEQUENCE:
#  wfe_new = defocus(wfe_orig, disp, fratio, lambda)
#INPUT PARAMETERS:
#  wfe_orig -- 2D image of WFE, in waves, associated with the optical system.
#  disp -- displacement (change in distance from aperture to detector), in 
#     same units as lambda.
#  fratio -- f-ratio of the optical system: effective focal length divided
#     by aperture diameter, where the aperture diameter equals the width of
#     the WFE image.
#  lam -- wavelength of light, in same units as distance. If included,
#     then it is assumed that the wfe input and output will be in waves.
#     If it is not given, then lambda is set to 1 (equivalently, we may
#     say that the WFE is assumed to be in the same distance units as
#     disp, or that the displacement is given in waves).
#OUTPUT PARAMETERS:
#  wfe_new -- wavefront error with a paraboloid of revolution added, to
#     simulate the desired defocus (specified by distance).
#OPTIONAL OUTPUT KEYWORDS:
#  wfe_change -- change in WFE (wfe_new - wfe_orig).
#  rmswfe -- RMS wavefront error.
#ALGORITHM:
#  We assume a large f-ratio, which results in a paraboloidal change in 
#  wavefront error given by:
#     WFE_change =  -0.5*(disp/lambda) * (x^2 + y^2) / fratio^2,
#  where x and y are coordinates centered in the aperture and normalized
#  to its width. If fratio < 4, a warning is thrown.
#BUGS:
#  There may be some ambiguity as to what is meant by a positive wavefront
#  error. If the convention of this program differs from the convention
#  that generated wfe_orig, then the sign of disp would effectively be
#  reversed. Might want to include a keyword for that. By the way, this
#  version has been tested to agree with an alternate calculation in
#  moses_defocus.
#MODIFICATION HISTORY:
#  2009-Jun-10  C. Kankelborg
#  2009-Jun-11  CCK reversed sign of WFE change based on testing with
#      moses_defocus.
#  2009-Aug-10  CCK Remove piston from wfe. Added rmswfe keyword.
#  2015-Oct-30  JP Modified piston removal line. See added comment.
#  2019 Jun 20  JTE translated into Python 3


import numpy as np

def defocus(wfe_orig,
            disp,
            f_ratio,
            lam = 1.0):
    
    if f_ratio < 4.0:
        print('WARNING: The approximation may be poor for fratio < 4.')
    
    #Create normalized aperture coordinates
    
    Nx, Ny = wfe_orig.shape()
    
    xlin = np.divide(np.arange(Nx),
                     Nx - 1.0) - 0.5
    
    #ylin = np.divide(np.arange(Ny) - 0.5 * (Ny - 1), Nx - 1.0)
    
    x = xlin 
    
    y = np.ones(Nx) # ylin
    
    #Calculate WFE
    wfe_change = -0.5 * np.divide(disp,
                                  lam) * np.divide(np.square([x, y]).sum(),
                                                   np.square(f_ratio))
    
    wfe = wfe_orig + wfe_change
    
    #remove piston, which is physically meaningless. Changed to only subtract the mean of the real part so that the imaginary part (the aperture) remains unchanged.
    wfe -= wfe.real.mean()
    
    rmswfe = np.sqrt(np.square(wfe).mean())
    
    return wfe, #rmswfe, wfe_change
    
