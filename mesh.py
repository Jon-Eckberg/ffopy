#NAME:
#  MESH
#PURPOSE:
#  Form a precise, scaled image of transmission from a standard Luxel filter mesh,
#  which is 70 lpi and approximately 82% transmission (per 
#  http://luxel.com/applications/xray-euv-filters/, accessed 2014-Jul-10).
#CALLING SEQUENCE:
#  image = mesh(inchsize)
#OUTPUT:
#  image = a square float image, with size given by inchsize (default 1").
#INPUT PARAMETERS:
#  INCHSIZE = size of the square filter image, in inches. Default = 1.0. 
#     If POWEROF2 is not set, INCHSIZE will be rounded to the nearest 1/70
#     inch, and that value will be returned. If POWEROF2 is set, then
#     the number of pixels in the image array is rounded to the next higher
#     power of 2, and INCHSIZE is adjusted accordingly. 
#OPTIONAL KEYWORD INPUTS:
#  POWEROF2: If set, then the mesh image size will be a power of 2 
#     (whatever is the next size up, based on inchsize).
#  IWFE: Normally, mesh puts ones in the openings and zeros on the gridlines.
#     If IWFE is set, then mesh takes the form of the IMAGINARY PART of the 
#     wavefront error, with zeros on the openings and 1e6 on the gridlines.
#OPTIONAL KEYWORD OUTPUTS:
#  PIXEL_SIZE_INCH = pixel size in inches.
#  TRANSMISSION = actual transmission value.
#METHOD:
#  The image is made from 21x21 cells. The mesh openings are 19x19 pixel 
#  squares, set to transmission 1.0, leaving a 2-pixel-wide grid with
#  transmission 0.0. Note that (19/21)^2 = 0.8186, which I take to be close
#  enough to the nominal 82%. This is the smallest whole number ratio that
#  yields the nominal transmission to within less than 1%. Note that the
#  result, image, always has pixel size 1/21/70 = 0.00068027212 inch.
#HISTORY:
#  2014-Jul-10  C. Kankelborg
#  2014-Jul-14  CCK fixed bug where output could have a non-integer number
#     of mesh cells. Also fixed errors in code comments.
#  2018-Apr-03  CCK added /powerof2 option.
#  2018-Apr-04  CCK added /iwfe option.
#  2019 Jun 20  JTE translated into Python 3

import numpy as np

def mesh(inchsize = 1.0,
         powerof2 = False,
         iwfe = False):

    #Geometrical setup. See notes under METHOD, above.
    pixels_per_cell = 21
    
    pixels_filled_per_cell = 19
    
    #Default filter size in inches. Modified later to actual.
    #Note that keyword_set argument need not be a keyword!
    
    lpi = 70.0 #lines (cells) per inch, fixed per Luxel spec
    
    if powerof2:
        
       Npixels = np.round(inchsize * lpi * pixels_per_cell)
          #Array size based only on inchsize.
       Npixels = np.power(2, np.ceil(np.log2(Npixels)))
          #Round up to the next larger power of 2.
       inchsize = np.divide(Npixels,
                            lpi*pixels_per_cell) #Make sure the size is consistent
          #with having an power of 2 pixels on a side.
    else: #default behavior: integer number of cells
       Ncells = np.round(inchsize * lpi) #number of mesh cells in the image.
          #We ensure that this is an integer by rounding here.
       Npixels = Ncells * pixels_per_cell #number of pixels on a side.
      
       inchsize = np.divide(Npixels,
                            lpi * pixels_per_cell) #Make sure the size is consistent
          #with having an integer number of cells.
    
    
    #Create 1D mesh.
    mesh1d = (np.arange(Npixels) % pixels_per_cell + 1) <= pixels_filled_per_cell
    
    #Create & return the 2D image.
    mesh2d = mesh1d # mesh1d
    
    if iwfe: #optional conversion to imaginary WFE
       mesh2d -= 1 #now we have -1 on the wires and 0 in the holes.
       mesh2d *= -1e6 #when squared, this will be 1e9 on the wires.
    
    #Calculate optional keyword outputs
    pixel_size_inch = np.power(lpi * pixels_per_cell, -1.0) #for optional output keyword.
    
    transmission = np.divide(mesh2d.sum(), mesh2d.size())
    
    #Return mesh image
    return mesh2d #,pixel_size_inch, transmission