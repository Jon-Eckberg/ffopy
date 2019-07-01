#NAME:
#   WFE2FRINGES
#PURPOSE:
#   Convert a map of wavefront error to an image of horizontal or vertical
#   fringes, as they would appear in an interferometer. It may be useless,
#   but it is cute.
#CALLING SEQUENCE:
#   fringes = wfe2fringes(wfe [, nfringes=nfringes] [, /vertical])
#INPUT PARAMETERS:
#   WFE = An Nx x Ny wavefront error map (in units of waves).
#OPTIONAL KEYWORD INPUTS:
#   NFRINGES = Number of fringes over the mirror surface. Default=4.
#   VERTICAL = If set, then make the fringes vertical (by default, they
#       are horizontal).
#OUTPUTS:
#   FRINGES = An image of the fringes you would see corresponding to the input WFE.
#MODIFICATION HISTORUY:
#   2008-OCT-30  C. Kankelborg
#   2008-Nov-07  CCK Made compatible with complex WFE maps (see psd2wfe.pro)
#   2019 Jun 21  JTE Translated into Python 3


import numpy as np

def wfe2fringes(wfe,
                nfringes = 4,
                vertical = False):
    
    Nx, Ny = wfe.shape
    
    #Construct perfect fringes, either vertical or horizontal.
    
    if vertical: #vertical fringes 
        fringes = nfringes * np.arange(Ny,Nx)/( (Nx-1)*(Ny-1) )
        
        fringes = fringes.T
        
    else: #horizontal fringes
        
        fringes = nfringes * np.divide(np.arange(Nx,Ny),  (Nx-1)*(Ny-1) )

    fringes += wfe  #Add wavefront error onto the perfect fringes

    return np.exp(2.0j * np.pi * fringes)
    
