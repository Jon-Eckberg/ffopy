#+
#NAME:
#  wfe2slope
#PURPOSE:
#  From a WFE array, compute an ISO 10110 compliant map of slope error. Units
#  of slope will be the ratio of the units of keyword lambda to the units 
#  of keyword dx, so these must match to get radians. If wfe2slope is provided
#  neither keyword, and the wfe is in waves (which we consider standard), then
#  the slope will be in waves per pixel. BY DEFAULT THE WAVEFRONT SLOPE IS
#  RETURNED, NOT THE UNDERLYING REFLECTIVE SURFACE. IF YOU WANT THE LATTER,
#  USE /SURF.
#CALLING SEQUENCE:
#  result = wfe2slope(wfe, integration_length, filter_length $
#        [, dx=dx][, lambda=lambda][, /surf])
#OPTIONAL KEYWORD INPUTS:
#  dx = wfe pixel size (default 1.0). The run of the calculated slope will have
#     the same units as dx.
#  lambda = wavelength (default 1.0). The rise of the calculated slope will
#     have the same units as lambda.
#  surf = if set, then return slope of the (reflecting) surface, rather than
#     wfe slope. This just divides the result by two!
#MODIFICATION SEQUENCE:
#  2015-Feb-21  C. Kankelborg
#  2015-Mar-09  CCK added /surf option.
#  2015-Mar-11  CCK corrected factor of two error in differential operators.
#-
def wfe2slope(wfe,
              integration_length,
              filter_length,
              dx = 1.0,
              lam = 1.0,
              surf=surf):

    
    Nx, Ny = wfe.shape()
    
    #CREATE DIFFERENTIATION OPERATORS (CONVOLUTION KERNELS)
    Npix = 1 + 2 * np.floor(0.5 * (1.0 + integration_length / dx))
       #Size of the convolution kernel we'll use for slopes.
       #This is the nearest odd integer to 1+integration_length/dx.
    
    integration_length_used = (Npix - 1) * dx
       #The integration length we are really using by forming
       #the convolution kernel with Npix pixels.
    
    d_dx = np.ones(Npix) #kernel for differentiation with respect to x.
    
    d_dx[0] = -1.0
    
    #d_dx[Npix - 1] = 1.0
    
    d_dx *= np.divide(lam, integration_length_used)
    
    d_dy = np.ones((1, Npix))  #kernel for differentiation wrt y.
    
    d_dy[0] = -1.0
    
    #CREATE PROPERLY PADDED, MASKED AND FILTERED WFE ARRAY
    Npad = Npix/2 + 1 #Leave room for the kernels to hang off the edge.
    
    wfe_big_real = np.full(Nx+Npad, Ny+Npad, np.nan) 
    
       #Real part. Pad with NaN
    wfe_big_real[0:Nx-1, 0:Ny-1] = wfe.real
       #Note that wfe_big has to be hit with convol(/edge_wrap) to take
       #advantage of padding on left & lower edges by wrapping.
    wfe_big_imaginary = np.full((Nx+Npad, Ny+Npad), 1e10)
    
    #This puts the padding outside the aperture, E field --> 0.
    wfe_big_imaginary[0:Nx-1, 0:Ny-1] = wfe.imag
    
    ss = np.where(wfe_big_imaginary > 1.0) 
       #Out-of-aperture indices, a somewhat arbitrary definition.
       #Where the imaginary part equals 0, the intensity is E^2 = 1.
       #Where the imaginary part equals 1, the intensity is E^2 = 3.5e-6.
    wfe_big_real[ss] = np.nan #Put NaN in out-of_aperture locations
    #FILTER THE REAL PART
    
    Nfilter = np.divide(filter_length, dx).round()
    
    wfe_big_real = smooth(wfe_big_real, Nfilter, /edge_wrap, /NaN)
       #This will result in some cross-talk between horizontal and vertical gradients
       #around the edge of a non-square mirror. However, this is the sort of filtering
       #that appears to be ISO standard.
       
    wfe_big_real[Nx:*,*] = np.nan #Restore the NaN padding that was smeared over by smooth.
    wfe_big_real[*,Ny:*] = np.nan #Restore the NaN padding that was smeared over by smooth.
    
    #COMPUTE GRADIENT
    slope_x = convol(wfe_big_real, d_dx, /edge_wrap)
    slope_y = convol(wfe_big_real, d_dy, /edge_wrap)
    
    if surf: 
        
       slope_x /= 2.0
       slope_y /= 2.0

    
    #COMPUTE RMS values
    #?? It would be even better to use imaginary part of wfe to weight the means ??
    slope_x_rms = np.sqrt(np.nanmean(np.square(slope_x)))
    slope_y_rms = np.sqrt(np.nanmean(np.square(slope_y)))
    slope_rms = np.sqrt(np.square([slope_x_rms, slope_y_rms]).sum())
       #I am assuming that there are similar numbers of
       #finite results in slope_x and slope_y, even though the
       #domains may not perfectly overlap!
    
    return slope_rms