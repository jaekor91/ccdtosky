# Configuration file for the program. All variables contained here must be
# specified by the user.

# File address for the ccd summary file
faddress = "../data/ccd/ccds-annotated-decals.fits" 

# Output directory
out_directory = "./"

# Parallelization parameters for computing statistics
num_cores = 1 # One by default. Non-negative integer. As an example, if you requested 8 cores, 
            # then you should set this number to8.

# Search radius range for each HEALPix center.
sepdeg = 0.336/2. 

# HEALPix parameters
Nside = 2**4 # Recommend 2**11 for accurate computation. 
			# WARNING: If more than 2**11, then compute time might be excessively long.
            # If less than 2**9, the approximation scheme used may not work as well.
NESTED = True # Use nested HEALPix division by default for histogramming.  
			# WARNING: Although the program should run correctly when False, do not change.

# This table is provided for your convenience.
# Nside 		# Pix  	      deg^2 	   arcsec^2   A_p/A_ccd    Targ # per pix
#  5 -         12,288     3.3571746      43,508,982     74.6039    8,057.219
#  6 -         49,152     0.8392936      10,877,245     18.6510    2,014.305
#  7 -        196,608     0.2098234       2,719,311      4.6627      503.576
#  8 -        786,432     0.0524559         679,827      1.1657      125.894
#  9 -      3,145,728     0.0131140         169,956      0.2914       31.474
# 10 -     12,582,912     0.0032785          42,489      0.0729        7.868
# 11 -     50,331,648     0.0008196          10,622      0.0182        1.967
# 12 -    201,326,592     0.0002049           2,655      0.0046        0.492
# 13 -    805,306,368     0.0000512             663      0.0011        0.123
# 14 -  3,221,225,472     0.0000128             165      0.0003        0.031
# Use Nside = 2^11. Pixel size is about 1 percent of the ccd size.


# Quantities of interests: For any quantity that the user is interested in, 
# user must specify a list of tuples (one per quantity) as in the following example.
import numpy as np
def func_nonsense(quantities, weights):
    """
    Dummy function for testing.
    """
    return np.sum(quantities*weights)

templates = [("Nexp","none", "sum"),
             ("airmass","galdepth_ivar", "mean"),
             ("ebv","galdepth_ivar", "mean"),
             ("seeing","galdepth_ivar", "mean"),
             ("avsky","galdepth_ivar", "mean"), # Is ADU to flux conversion necessary?
             ("airmass","galdepth_ivar", func_nonsense)]

# More generally,
#
# ("ccd_quantity", "weight", operation)
#
# - "ccd_quantity": str. For examlpe, "seeing", "airmass", "ebv", etc.   
# - "weight": str. There are two choices. "none" does not use any kind of weighting scheme. 
#	If "galdepth_ivar", then galactic depth inverse variance in invesre nanomaggies square 
#	is used as weights.   
# - operation: str or Python function. If not one of the operations below, then operation 
#	must be a Python function (e.g., func_nonsense) that takes exactly one (if "weight" = "none") 
#	or two (if "weight" = "galdepth_ivar") arguments. The function must be defined within the 
#	config file.
#     - "mean": If "weight" = "galdepth_ivar", then weighted average is computed.
#     - "min"
#     - "max"
#     - "median"
#     - "sum"
#     - "std": numpy.std
#  
# Note that the output quantity is computed for each band seperately ("g", "r" and "z"). 
# The output recarray will have 3*n+3 columns, where n is the number of quantities specified. 
# The array contains HEALPix pixels and the corresponding ra/dec's.







