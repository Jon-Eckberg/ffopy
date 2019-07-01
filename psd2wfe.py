#NAME:
#   PSD2WFE
#PURPOSE:
#   Convert a power spectral density for polishing errors in a
#   reflective optic to a randomly generated wavefront error (WFE) map.
#   The WFE is generated from the power spectrum using random phases
#   and amplitudes. It is *not* periodic. Note from the documentation
#   below that the PSD runs from zero to the Nyquist frequency in units
#   of exactly *half* the fundamental. Note that the output WFE map is
#   complex (see keyword CIRCULAR). The imaginary part could be used
#   to represent apodization or nonuniform illumination of the aperture
#   if you like.
#CALLING SEQUENCE:
#   wfe = psd2wfe(psd [, rms=rms | pv=pv])

#INPUT PARAMETERS:
#   PSD = An N-element power spectral density for the polishing errors.
#       The period associated with the ith element of PSD is 2(N-1)/i
#       pixels in the ouptut N x N wafefront error map. One consequence
#       of this is that there is noplace to specify power in the frequency 
#       interval between the 1D Nyquist frequency and the diagonal 
#       Nyquist frequency. The PSD will be assumed flat in this range.

#OPTIONAL KEYWORD INPUTS:
#   RMS = The (scalar) RMS wavefront error. This is used to scale the
#       result, which eliminates any ambiguity about the normalization
#       of the Fourier transform (and hence of the PSD).
#   PV = Alternative to RMS, may specify peak-to-valley.
#   CIRCULAR = If set, then a circular aperture is inscribed in the square
#       WFE array. The region outside the circle is assigned the value
#       complex(0, 1e9). This results in a zero E-field later, which is
#       appropriate for regions outside the aperture. 
#   VERBOSE =  If set, then print out some informative statistics.
#   OVERSIZE = factor by which to oversize the WFE map, by comparison to the
#       diameter. Default=1.

#OUTPUTS:
#   WFE = a randomly generated, N x N wavefront error map that is consistent
#         with the input PSD and RMS or PV. If neither RMS nor PV are specified,
#         then the WFE will have a PV of 0.25 (a notional quarter-wave mirror).

#MODIFICATION HISTORY:
#   2008-OCT-29  C. Kankelborg
#   2008-Nov-03  CCK Incorporated removal of best-fit spherical wave
#       before normalizing the WFE.
#   2008-Nov-07  CCK Added keyword option /CIRCULAR, and thereby 
#       introduced the concept of a complex WFE map. The imaginary
#       part represents attenuation of the electric field amplitude.
#   2008-Nov-10  CCK fixed bug--- np.sqrt(PSD) to get amplitudes. Added 
#       diagnostic printing of Strehl ratio. Added VERBOSE keyword
#       so that diagnostic messages are suppressed by default.
#  2009-Oct-16  CCK added OVERSIZE keyword to allow room for the beam to be
#       stretched, e.g. by oblique reflection from a grating.
#  2018-Apr-24  CCK added KEEP_DEFOCUS keyword to suppress the use of 
#        FIT_SPHERICAL_WAVE, which is time consuming. This speeds things up when
#        the routine is being called with small or zero RMS, perhaps just to create 
#        an aperture map.
#  2019-Jun-19 JTE translated function to python 3

import numpy as np

def psd2wfe(psd,
            rms = 0,
            pv = 0.25, #quarter wave nominal
            circular = False,
            oversize = 1.0,
            keep_defocus = False,
            verbose = False):
  
   
    N = psd.size
    
    frequencies = np.fft.fftfreq(2 * N) < N  #2n_psd x 2n_psd array of frequencies
    
    ft_power = interpolate(psd, frequencies, cubic= - 0.5) > 0
    
     #threshold is to eliminate negatives that creep in from the interpolation.
    ft_amplitudes = np.sqrt(ft_power)
    
    #amplitude is the square root of power.
    ft_phases = 2.0 * np.pi * np.random.rand(2 * N, 2 * N)
    
    ft_complex = ft_amplitudes * np.exp(np.sqrt(1.0j) * ft_phases)
    
    wfe_big = np.fft.fft2d(ft_complex).real
    
    wfe1 = wfe_big[:N-1, :N-1] #Cropping eliminates the periodic boundary conditions.
    
    #Now subtract best-fit spherical wave
    if keep_defocus:
        
        spherewave = 0.0
        
    else:
        
        spherewave = fit_spherical_wave(wfe1)
        
    wfe = np.complex(wfe1 - spherewave)
    
#   Henceforth WFE will be complex. This has many uses, including
#   the representation of apetures, apodization, and nonuniform
#   illumination.
    if verbose:
        
        if not keep_defocus:
        
            print('Subtracted spherical component PV = {}'.format(spherewave.ptp()))
    
        print('Residual WFE before normalization is PV = {}'.format(wfe.ptp()))
        
    #Introduce aperture cropping, if any
    if circular:
       
        #fix here
    
        axis = np.linspace(-N//2, N//2, N)
        
        radius = np.sqrt(axis**2 + axis[:, np.newaxis]**2)
        
        aperture = np.where(radius < (oversize * 0.5 * N))
        
        crop = np.where(radius >= (oversize * 0.5 * N))
       
        wfe[crop] = 1e9j
    
    else:
        
        aperture = np.mgrid[:N, :N]
     
    #Normalize wfe (any apodization or nonuniform aperture features should
    #be inserted AFTER this step, since we would not want them to interfere
    #with the scaling of the mirror surface figure)
    
    wfe[aperture] -= wfe[aperture].mean() #set to mean of zero.

    if rms:
        
        wfe[aperture] *= np.divide(rms, wfe[aperture].std()) 
          #scale to desired RMS wavefront error

    else:
        
        wfe_range = wfe[aperture].ptp()
        
        wfe[aperture] *= np.divide(pv,
                                   wfe_range)
    
    if verbose:
        
        print('WFE realized  PV = {}'.format(wfe[aperture].ptp()))
        
        print('WFE realized RMS = {}'.format(wfe[aperture].std()))
        
        print('Strehl = {}'.format(wfe2strehl(wfe)))
        
    return wfe