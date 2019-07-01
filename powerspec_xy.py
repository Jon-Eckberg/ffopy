#NAME:
#  POWERSPEC_XY
#PURPOSE:
#  Estimate horizontal and vertical power spectra of an image. The horizontal power
#  spectrum, for example, is calculated as an average of the power spectra of the rows.
#  Moreover, the power spectrum of each row is estimated by averaging over a sequence
#  of short segments along the row. Columns are treated similarly. The normalization
#  convention is the same as with periodogram_fft.pro.
#CALLING SEQUENCE:
#  powerspec_xy, image, psx, psy [, k=k] [, nfreq=nfreq]
#INPUT PARAMETERS:
#  image = the image for which the power spectra are to be calculated.
#OUTPUT PARAMETERS:
#  psx,psy = the horizontal and vertical power spectra, respectively.
#OPTIONAL KEYWORD INPUTS:
#  nfreq = Number of elements in each segment, and also the number of
#     elements in each of the output power spectra. Default is 64. 
#     The smaller nfreq, the lower the spectral resolution but the 
#     better the noise suppression.
#  dx = sample interval. Default=1. If included, then k will be converted into units
#     of cycles (or radians, if that keyword is set) per unit distance (or time),
#     corresponding to the units of dx.
#  radians = if set, then let k be in radians (rather than cycles) per whatever.
#  retain_dc = if set, then retain the DC offset when calculating the power spectrum.
#     Ordinarily, the DC offset is subtracted before calculating the spectrum, because
#     the Hanning window tends to make the offset bleed into the nearest neighbors of
#     pixel [0,0] of the 2D power spectrum.
#OPTIONAL KEYWORD OUTPUTS:
#  k = Array of wavenumbers (by default, in cycles per pixel).
#DEPENDENCIES:
#  powerspec.pro
#MODIFICATION HISTORY:
#  2013-Sep-05  C. Kankelborg
#  2013-Oct-03  CCK Provided checking for finite power spectra, so that infinities
#     can now be excluded from the average power spectrum.
#  2013-Oct-03 bug report:
#IDL> .run test_random_image
#% Compiled module: $MAIN$.
#% Compiled module: RANDOM_IMAGE.
#% Compiled module: K_ARR2D.
#% Compiled module: K_ARR.
#% Compiled module: PERIODOGRAM_FFT.
#% Compiled module: HANNING.
#% Compiled module: POWERSPEC2D.
#% Program caused arithmetic error: Floating divide by 0
#IDL> window, 10
#IDL> data = image - mean(image)
#IDL> powerspec_xy, image, psx, psy
#% Compiled module: POWERSPEC_XY.
#% Compiled module: POWERSPEC.
#IDL> print, mean(psx)
#      4185.73
#IDL> print, mean(psy)
#      8.22428
#IDL> print, mean(data^2)
#      4418.29
#******************************************************************
#
#powerspec_xy, image, psx, psy, k=k, nfreq=nfreq, $
#   dx=dx, radians=radians, retain_dc=retain_dc

Nx, Ny = image.shape()

if (nfreq.size() != 1):
    nfreq = 64

psx = fltarr(Nfreq)
Nvalid = Ny #A priori, assume we'll have Ny valid power spectra.
for y in range(Ny):
    
   ps = powerspec(reform(image[:, y]),
                  nfreq,
                  dx = dx,
                  kx = k,
                  radians = radians,
                  retain_dc = retain_dc)
   
   if finite(ps[0]):
      psx += ps
   else:
      Nvalid -= 1 #one less valid power spectrum than we thought!
   
psx /= Nvalid #take average over good power spectra.

psy = fltarr(Nfreq)

Nvalid = Nx #A priori, assume we'll have Nx valid power spectra.

for x in range(Nx):
    
   ps = powerspec(reform(image[x]),
                  nfreq,
                  dx=dx,
                  kx=k,
                  radians=radians,
                  retain_dc = retain_dc)
   
   if finite(ps[0]):
      psy += ps
   else:
      Nvalid -= 1 #one less valid power spectrum than we thought!


psy /= Nvalid #take average over good power spectra.
