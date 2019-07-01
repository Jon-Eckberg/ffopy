#NAME:
#  periodogram_fft
#PURPOSE:
#  Perform the simplest reasonable power spectral estimate, using the FFT
#  and a Hann window. The normalization is such that the mean value of the
#  result is approximately the mean squared value of the data. Works 
#  for 1D or 2D inputs. It should be noted that the expected error of this
#  estimate of the power spectrum is about 100% (see Numerical Recipes).
#  However, this routine is needed by the much better estimators powerspec
#  and powerspec2d.
#CALLING SEQUENCE:
#  result = periodogram_fft(data)
#INPUTS:
#  data = 1d or 2d array for which the power spectrum is desired. It is not assumed
#     that the data is periodic# therefore, the data will be Hann windowed. It
#     is also not assumed that the image is real. That's why the normalization
#     is such that the sum must be performed over both positive and negative k.
#OPTIONAL KEYWORD INPUTS:
#  retain_dc = if set, then retain the DC offset when calculating the FFT.
#     Ordinarily, the DC offset is subtracted before calculating the spectrum, because
#     the Hanning window tends to make the offset bleed into the nearest neighbors of
#     the DC element of the FFT.
#NOTES:
#  Use k_arr or k_arr2d to get the k-values for each element in the result.
#HISTORY:
#  2013-Aug-20  CCK
#  2013-Sep-06  CCK fixed bug in 1D case (was not using data_corrected).
#  2013-Oct-03  CCK fixed so that 1xN arrays will be treated as 1D.
#  2019 Jun 25  JTE translated into Python 3

import numpy as np

def periodogram_fft(data,
                    retain_dc = True):

    if retain_dc:
    
        data_corrected = np.squeeze(data)
    
    else:
    
        data_corrected = np.squeeze(data - data.mean())
    
    #The corrected data has had trailing zero dimensions eliminated. This way,
    #for example, 1xN arrays will be trated as 1D rather than 2D.
        
    Ndim = data_corrected.ndim
    
    Nx, Ny = data_corrected.shape
    
    if (Ndim > 2):
           print('Data must be a 1d or 2d array.')
    
    if (Ndim == 2): #It's an image.
        
        power = np.square(np.abs(np.fft.fft2(data_corrected * np.hanning(Nx) * np.hanning(Ny)[:, np.newaxis]))) * 7.1111111 * Nx * Ny
          #The factor 7.1111111 is one over the mean square of the Hann window.
          
    else: #It's a time series or other 1d array (or else it's an error!).
        
       power = np.square(np.abs(np.fft.fft2(data_corrected * np.hanning(Nx)))) * 2.6666667 * Nx
          #The factor 2.6666667 is one over the mean square of the Hann window.
    
    return power