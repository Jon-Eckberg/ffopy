#NAME:
#  wfe2strehl
#PURPOSE:
#  Convert a wavefront error map to a Strehl ratio. This is defined
#  as the fraction of diffraction-limited intensity at the
#  diffraction-limited (gaussian) focus. For the result to be correct,
#  the best-fit sphere must have been subtracted from the WFE. All my
#  functions that generate WFE maps should do this automatically 
#  (e.g. psd2wfe).
#CALLING SEQUENCE:
#  strehl = wfe2strehl(wfe)
#INPUT PARAMETERS:
#  wfe = wavefront error map in waves.
#OUTPUT PARAMETERS:
#  strehl = the strehl ratio. According to the Marechal criterion,
#     strehl > 0.8 is diffraction limited.
#ALGORITHM:
#  I use the exact form of the Strehl ratio given by Wyant ch.1, 
#  eq. 64. See:
#     http://www.optics.arizona.edu/jcwyant/zernikes/Zernikes.pdf
#MODIFICATION HISTORY:
#  2008-Nov-09  C. Kankelborg

import numpy as np

def wfe2strehl(wfe):
    
    Efield = np.exp(2.0j * np.pi * wfe)
    
    strehl = np.square(np.divide(np.abs(Efield.sum()),
                                 np.abs(Efield).sum()))
    
    return strehl