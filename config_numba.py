# Configuration file for the program. All variables contained here must be
# specified by the user.

# File address for the ccd summary file
faddress = "../data/ccd/ccds-annotated-decals.fits" 

# Output directory
out_directory = "./outputs/decals_Nside9/"

# If True, use numba optimized function in step 6. Otherwise, use 
# scipy.stats.binned_statistic
use_numba = True

# Numba number of threads. Irrelevant if use_numba = False. 
num_thread = "2"
# Must be a str. Python binned statistic function is defined with numba jit 
# decorator with nogil=True, meaning that the function may gain performance
# with a greater number of threads.

# Use spherematchDJE code to perform spherematching. A lot faster than astropy
# routine.
DJE = True
# WARNING: The C-program assumes that the user's machine treats int as 32-bits
# integer and double as 64-bits floating point number. If this is not the case
# use DJE = False. 

# Search radius range for each HEALPix center.
# WARNING: Ideall, this number should not be changed.
sepdeg = 0.336/2. 

# HEALPix parameters
Nside = 2**7    # Recommend 2**11 for accurate computation. 
			# WARNING: If more than 2**11, then compute time might be excessively long.
            # If less than 2**9, the approximation scheme used may not work as well.
NESTED = False # Use nested HEALPix division by default for histogramming.  

# This table is provided for your convenience.
# Nside 		# Pix  	deg^2 	   arcsec^2   A_p/A_ccd    Targ # per pix
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
templates = [
             ("Nexp","none", "sum"),
             ("airmass","galdepth_ivar", "max"),
             ("airmass","galdepth_ivar", "min"),
             ("airmass","galdepth_ivar", "mean"),
             ("ebv","galdepth_ivar", "max"),
             ("ebv","galdepth_ivar", "min"),
             ("ebv","galdepth_ivar", "mean"),             
             ("seeing",	"galdepth_ivar", "max"),
             ("seeing", "galdepth_ivar", "min"),
             ("seeing", "galdepth_ivar", "mean"),             
             ("galdepth","none", "none")
             ]

# More generally,
#
# ("ccd_quantity", "weight", operation)
#
# - "ccd_quantity": str. For examlpe, "seeing", "airmass", "ebv", etc. If "Nexp" is used
#					then unweighted sum is forced.
# - "weight": str. There are two choices. "none" does not use any kind of weighting scheme. 
#			If "galdepth_ivar", then galactic depth inverse variance in invesre nanomaggies square 
#			is used as weights.   
# - operation: str. The following statistics are available:
#     - "mean": Compute the (weighted) mean of values for points within each bin. Empty bins will be 
#				represented by NaN. 
# 	  - "median": compute the median of values for points within each bin.
#		 		Empty bins will be represented by NaN. "weight" is ignored.
#                       Only available for scipy option for the moment.
# 	  - "min": compute the minimum of values for points within each bin.
#				 Empty bins will be represented by NaN. "weight" is ignored.
# 	  - "max": compute the maximum of values for point within each bin. 
#				Empty bins will be represented by NaN. "weight" is ignored.
# 	  - function: a user-defined function which takes a 1D array of values, 
#				and outputs a single numerical statistic. This function will 
#				be called on the values in each bin. Empty bins will be represented 
#				by function([]), or NaN if this returns an error.
#  
# Note that the output quantity is computed for each band seperately ("g", "r" and "z"). 
# The output recarray will have 3*n+3 columns, where n is the number of quantities specified. 
# The array contains HEALPix pixels and the corresponding ra/dec's.