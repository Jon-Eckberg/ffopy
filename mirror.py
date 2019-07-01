#NAME:
#  mirror
#PURPOSE:
#  Use physical parameters to generate PSD and WFE for a mirror.
#  The polishing process is modeled by a power law PSD from 1/2
#  period per diameter down to a cutoff.
#CALLING SEQUENCE:
#  wfe = mirror(diameter, rmswfe [, N=N] [, hole=hole] 
#     [, pindex=pindex] [, pcutoff=pcutoff] 
#     [, psd=psd] [, period=period]
#INPUT PARAMETERS:
#  diameter = diameter of the beam on the mirror (not the mirror
#     itself) in physical units (usually mm). Optionally, diameter
#     may be a 2-element array giving the width and height of an ellipse.
#OPTIONAL INPUT PARAMETERS:
#  rmswfe = RMS wavefront error (floating point number). If not supplied,
#     or if zero is given, then the default behavior of psd2wfe will be 
#     used (1/4 wave PV).
#OUTPUT PARAMETERS:
#  wfe = wavefront error (WFE) for a particular, randomized
#     realization with the desired PSD.
#OPTIONAL KEYWORD INPUTS:
#  N = number of points (integer). The resulting PSD will have N
#     elements, and the WFE map will be NxN. Default=256.
#  circular = if set, use a circular aperture (passed straight
#     through to psd2wfe).
#  pindex = power-law index of the power spectral density (PSD)
#     of the polishing errors. PSD goes as 
#           PSD = k**-pindex. 
#     Typical values are 2.1-3.0. Default = 2.5.
#  pcutoff = The shortest period that will have nonzero power.
#     Default = 4 (e.g. 4mm if diameter is in mm).
#  hole = Diameter of concentric hole to be drilled in mirror. Uses same
#     units as diameter and pcutoff.
#  renormalize = If set, then renormalize to get the specified rmswfe over
#     the clear aperture, excluding the "hole" if specified. If not set,
#     then rmswfe will pertain to everything within a circle with diameter
#     equal to the maximum value of the diameter argument.
#  oversize = factor by which to oversize the WFE map, by comparison to the
#     *minimum* value in the diameter array. I scale from the minimum because
#     the minimum beam dimension is more likely to remain constant as
#     you go through a multi-element system like a spectrograph. Default=1.
#  keep_defocus = if set, then do not subtract the best-fit spherical wave. This
#     is just a pass-thru to psd2wfe.
#OPTIONAL KEYWORD OUTPUTS:
#  psd = The N-element PSD array. Note that the scale of the PSD
#     values is arbitrary.
#  period = N-element array of period values, in the same units
#     as diameter and pcutoff.
#  aperture = index array of WFE elements within the clear aperture.
#     Especially handy for annular mirrors (see CIRCULAR & HOLE keywords).
#  dx = pixel size in same units as the diameter.
#  width = width of the wfe map in units of the mirror diameter.
#MODIFICATION HISTORY:
#  2008-Nov-09  C. Kankelborg
#  2008-Nov-10  CCK Added HOLE and APERTURE keywords. Modified default 
#       pindex to 2.5 based on SAO input. Improved documentation.
#  2008-Nov-17  CCK added RENORMALIZE keyword to correct rmswfe when HOLE
#       is specified.
#  2009-Oct-16  CCK added  OVERSIZE keyword to allow room for the beam to be
#       stretched, e.g. by oblique reflection from a grating. Also added
#       WIDTH output keyword. Also modified to allow an elliptical mirror
#       by specifying a 2-element diameter.
#  2015-Feb-06  CCK corrected minor error in doc header PURPOSE section.
#  2015-Feb-20  CCK fixed documentation of dx units (same units as diameter,
#        rather than units "of" the diameter). Then discovered that dx had
#        never actually been implemented, so I fixed that!
#  2015-Mar-11  CCK clarified documentation for default behavior if rmswfe=0.
#  2018-Apr-24  CCK added keep_defocus keyword (see psd2wfe.pro).
#  2019-Jun-19  JTE translated to Python 3

import numpy as np

def mirror(diameter,
           rmswfe,
           N = 256,
           pindex = 2.5,
           pcutoff = 4.,
           psd = None,
           circular = False,
           period = period,
           aperture = aperture,
           hole = 0.0,
           renormalize = False,
           oversize = 1.0,
           keep_defocus = False):

    #Default keyword values
    
    width = diameter.min() * oversize
    
    dx = np.divide(diameter,
                   N - 1) 
    #CCK 2015-Feb-20# dx keyword finally implemented!
    
    k = np.arange(N) #wavenumber scale, arbitrary units
    
    period = np.divide(2.0 * width,
                       k) #period in same units as diameter
    
    #NOTE:
    #In keeping with psd2wfe, an N-element PSD leads to an NxN WFE.
    #The following table may help to clarify the meaning of the 
    #frequency elements:
    #
    #  k        interpretation
    #  ----------------------------------
    #  0        DC.
    #  1        1/2 period per map width
    #  2        1 period per map width
    #  ...
    #  N-1      Nyquist frequency
    
    
    #Implement power law PSD
    psd = np.power(k, -pindex)
    
    psd[0] = 0
    
       #first k with nonzero power is 1 per diameter
    
    #Implement cutoff frequency
    ss = np.where(period < pcutoff)
    
    if (ss[0] != -1):
        
        psd[ss] = 0
    
    #Calculate and return WFE
    wfe = psd2wfe(psd,
                  rms = rmswfe,
                  circular = False,
                  oversize = False,
                  keep_defocus = False)
    
    #Drill hole if desired
    if hole > 0.0:
        
        radius = np.divide(N * hole, 2.0 * width).round() #calculate radius in pixels

        axis = np.linspace(-N/2, N/2, N)

        ########## fix dist        
        radii = np.sqrt(axis**2 + axis[:, np.newaxis]**2)
        
        ss = np.where(radii < radius) #find where the hole lives
        
        wfe[ss] = 1e9j   #drill the hole
    
    
    #Make elliptical aperture if desired
    if diameter.ndim == 2:
        
       xy = np.divide(width * (np.arange(N) - N//2), N) #x or y linear coordinate
       
       x = xy # replicate(1.0,N) #x coordinate array
       
       y = np.ones(N) # xy #y coordinate array
       
       ss = np.where(4.0 * np.square([x / diameter[0], y / diameter[1]]).sum() > 1.0 )
       
       wfe[ss] = 1e9j
    
    aperture = np.where(wfe.imag < 1.0 )
    
    if renormalize:
    
        wfe[aperture] *= np.divide(rmswfe, wfe[aperture].std())
        
    return wfe