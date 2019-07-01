#+
#NAME:
#  RMSWFE
#PURPOSE:
#  Given a complex map of the wavefront error, calculate the RMS wavefront error
#  with proper weighting of the aperture. NO ATTEMPT IS MADE TO SUBTRACT THE MEAN
#  WAVEFRONT ERROR, TIP, TILT, DEFOCUS, OR ANY OTHER TERMS.
#CALLING SEQUENCE:
#  result = rmswfe(wfe, strehl=strehl, approx_strehl=approx_strehl)
#OUTPUT:
#  The result is a floating point scalar RMS wavefront error in waves (assuming
#  wfe input is wavefront error in waves).
#INPUTS:
#  WFE = complex wavefront error. The real part is wavefront error in waves. The 
#     imaginary part is an amplitude. One calculates the scalar electric field as 
#     E = exp(i * 2 * !pi * WFE). Typically, we use an imaginary part of 0 in the
#     clear aperture of the optic, leading to abs(E) = 1. You might suppose that
#     a figure error input would result in an rms figure error output# that is correct
#     so long  as the imaginary part is encoded correctly for the aperture weighting.
#OPTIONAL OUTPUT KEYWORDS:
#  STREHL = The strehl ratio, calculated exactly from summing the field over the aperture.
#  APPROX_STREHL = Strehl ratio, calculated using the extended Marechal's approximation.
#HISTORY:
#  2013-May-01  C. Kankelborg
# 2019 Jun 21  JTE Translated into Python 3

import numpy as np

def rmswfe(wfe):
    
    E_real = np.exp(2.0j * np.pi * wfe).sum() #E-field at PSF center
    
    amplitude = np.exp(-2.0j * np.pi * wfe.imag) #E-field amplitude at aperture (with no OPD)
    
    E_ideal = amplitude.sum() #E-field at PSF center without aberrations

    strehl = np.power(np.divide(np.abs(E_real), np.abs(E_ideal)))
    
    rmswfe = np.sqrt(np.divide(amplitude * np.square(wfe.real).sum(), E_ideal))
    
    approx_strehl = np.exp(-np.square(2.0 * np.pi * rmswfe)) #extended Marechal's approximation

    return rmswfe#, strehl
