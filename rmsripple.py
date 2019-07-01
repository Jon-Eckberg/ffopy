#+
#NAME:
#  rmsripple
#PURPOSE:
#  Given a complex map of the wavefront error, calculate the RMS wavefront error
#  over a specified wavenumber band. This is useful for calculating RMS ripple,
#  i.e. the waviness over mid spatial frequencies. This routine relies on 
#  powerspec2d(), which uses Welch's method. Since the WFE image is broken up into
#  small pieces and the mean is subtracted from each, THIS ROUTINE IS NOT
#  APPROPRIATE FOR CALCULATING RMSWFE. For that task, use rmswfe(). This routine
#  is only appropriate for a band-limited calculation using a minimum frequency
#  significantly greater than the fundamental. This routine is not much good
#  for situations in which significant vignetting or apodization is represented 
#  in the WFE, because the imaginary part of the WFE is used only to mask the
#  clear aperture (aperture is considered open wherever imaginary(wfe) < 1).
#CALLING SEQUENCE:
#  result = rmsripple(wfe, ripple_period, ripple_bandwidth, dx=dx [, /surf])
#OUTPUT:
#  The result is a floating point scalar RMS wavefront error in the same
#  units as the wavelength keyword lambda. If lambda is not supplied,
#  the units will be the same as the input WFE map (nominally, waves).
#INPUTS:
#  WFE = complex wavefront error. The real part is wavefront error in waves. The 
#     imaginary part is an amplitude. One calculates the scalar electric field as 
#     E = exp(i * 2 * !pi * WFE). Typically, we use an imaginary part of 0 in the
#     clear aperture of the optic, leading to abs(E) = 1.
#  ripple_period = longest period for the band-limited measurement. This is
#     in the same units as dx. If dx is not supplied, the units will be pixels.
#     This value is also used to set the size of the sub-windows that will be
#     used for calculating the PSD in subroutine powerspec2d().
#  ripple_bandwidth = a floating point constant greater than one. This will be
#     the ratio of longest to shortest periods for the band-limiting filter.
#OPTIONAL KEYWORD INPUTS:
#  dx = pixel size, in whatever units are desired (mm is typical). Default=1,
#     which puts the units in pixels. The same units will be assumed for 
#     interpreting ripple_period.
#  surf = if set, then calculate the ripple of the (reflective) SURFACE.
#     This just divides the result (including PSD, if you output that) by two!
#  lambda = wavelength, default=1.
#OPTIONAL KEYWORD OUTPUTS:
#  PSD = Power spectral density of input WFE (or of the underlying reflective
#     surface, if /surf is set)
#  K = Array of wavenumbers, in cycles per (unit of dx).
#HISTORY:
#  2015-Feb-21  C. Kankelborg
#  2015-Mar-09  CCK  Added /surf option. Added optional outputs PSD, K.
#  2019 Jun 20  JTE translated into Python 3

import numpy as np

def rmsripple(wfe,
              ripple_period,
              ripple_bandwidth,
              k,
              dx = 1.0,
              lam = 1.0, #lambda
              surf = False):

    kmin = np.power(ripple_period, -1.0)
    
    if (ripple_bandwidth <= 1.0):
        
        print('ERROR: ripple_bandwidth must be >1.')
        
    kmax = ripple_bandwidth * kmin
    
    #Mask aperture
    wfe_real = wfe.real
    
    ss = np.where(wfe.imag > 1.0) 
       #Out-of-aperture indices, a somewhat arbitrary definition.
       #Where the imaginary part equals 0, the intensity is E^2 = 1.
       #Where the imaginary part equals 1, the intensity is E^2 = 3.5e-6.

    wfe_real[ss] = np.nan #Put NaN in out-of_aperture locations
    
    #Compute PSD
    Npix = np.ceil(np.divide(ripple_period, dx)) #size of subimage for powerspec2d.
       #We want the size of this subimage to be greater than or
       #equal to the ripple period in pixels.
    psd = powerspec2d(wfe_real,
                      Npix,
                      dx = dx,
                      k_magnitude = k)
    
    #Filter PSD by zeroing out-of-band elements.
    ss = np.where(k < kmin)
    
    if (ss[0] != -1):
        psd[ss] = 0.0
    
    ss = np.where(k > kmax)
    
    if (ss[0] != -1):
        psd[ss] = 0.0
        
    psd[0,0] = np.nan #ensure the DC pixel is not counted in the RMS calculation.
    
    if surf:
        psd /= 2.0
    
    rms_ripple = np.sqrt(np.nanmean(psd))
    
    rms_ripple *= lam #convert to the same units as lambda.
    
    return rms_ripple