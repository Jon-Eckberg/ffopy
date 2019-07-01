#NAME:
#   WFE2PSD
#PURPOSE:
#   Produce an estimate of the power spectral density (PSD) of a wavefront
#   error map (WFE).
#CALLING SEQUENCE:
#   psd_est = wfe2psd(wfe)
#ALGORITHM:
#   A periodogram approach is used with a Hanning window. Two estimates
#   are produced and averaged: one using a vertical chord, and one a
#   horizontal chord. Compare Walsh et al. 1999, Appl. Opt. Vol. 38, No. 22, 
#   p.4790.
#MODIFICATION HISTORY:
#   2008-Nov-10  C. Kankelborg
#   2019 Jun 21   JTE translated into Python 3

import numpy as np

def wfe2psd(wfe):
    
    wfe_size = wfe.shape()
    N = wfe_size[0]
    if wfe_size[1] != N:
        print('WFE not square. Exiting.')
        return
    
    window = np.hanning(N)
        #windowing function for estimate of PSD
    
    chord1 = wfe[N//2].real
    chord2 = wfe[:, N//2].real
    
    c1fft = np.fft.fft(chord1 * window)
    c2fft = np.fft.fft(chord2 * window)
    
    psd = np.square(np.abs(c1fft) + np.abs(c2fft))
    
    #IDEA: ADD FEATURE TO FIT A POWER LAW TO THE PSD
    
    return psd[:N//2]
