# Configuration file for the program. All variables contained here must be
# specified by the user.

# File address for the ccd summary file
faddress = "../data/ccd/ccds-annotated-decals.fits" 

# Output directory
out_directory = "./"

# This feature is current disabled.
# # Parallelization parameters for computing statistics
# num_cores = 1 # One by default. Non-negative integer. As an example, if you requested 8 cores, 
#             # then you should set this number to8.

# Search radius range for each HEALPix center.
# WARNING: Ideall, this number should not be changed.
sepdeg = 0.336/2. 

# HEALPix parameters
Nside = 2**9 # Recommend 2**11 for accurate computation. 
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

# Number of pixels in each unit before matching to ccd_centers. (Slows down
# the computation by ~2 but reduces memory requirement.)
Nside_kdtree = 2**7
# WARNING: Must be in powers of 2.
# A guide on the choice of Nside_kdtree:
# The expected number of matches between the ccd centers and pixel centers 
# can be estimated as follows: DESI covers about the third of the sky, so 
# the number of HEALPix pixels that overlap with ccd regions is about a 
# third of the total number.(For DR3, it's actually 18 percent.) There 
# will be num_ccd_avg number of ccds that corresponds to each pixel on 
# average. (For DR3, num_ccd_avg ~6.) Putting the two figures together 
# this means we expect about (num_pix * num_ccd_avg/3.) matches between ccds
# and the pixels. As an example, for DR3 with Nside = 2^11 case, we would 
# expect 54M = (50M * 6 * 0.18) matches. That's 432 MB of memory for numpy arrays
# for match indices. For the compelete DESI imaging survey, we would instead
# expect (as an upper bound) 330M = (50M * 20 * 1/3.) matches. That translates
# to almost 2.8 GB of memory. More importantly, the kd-tree algorithm that the program implements
# requires RAM memory that grows proportinate to the number of input pixels.
# For DR3, Nside = 2^11 case, the requirement is over 15GB! So to be conservative
# in memory requirement, we opt to divide the computation into several units 
# with each unit containing hp.nside2npix(Nside_kdtree) pixels. In Nside = 2^9 case, 
# kd-tree takes up about 3 GB (or less) and since the computation 
# is done on chunks, there is no worry of running out of RAM memory at least
# for computing pixel-ccd matches. However, this does mean the computation
# is done in 16 = (2^11/2^9) chunks, which will add to overhead. 
# If the number of matches grows by another factor of 10, say, because more 
# data was taken then that translates to 28 GB of memory just for indices!
# That would require saving the matched indices into files, out of RAM,
# which my current implementation does not do.

# Quantities of interests: For any quantity that the user is interested in, 
# user must specify a list of tuples (one per quantity) as in the following example.
templates = [
			 # ("Nexp","none", "sum"),
             # ("airmass","galdepth_ivar", "min"),
             # ("airmass","none", "mean"),             
             # ("airmass","galdepth_ivar", "mean"),
             # ("ebv","galdepth_ivar", "mean"),
             # ("seeing",	"galdepth_ivar", "mean"),
             ("avsky","galdepth_ivar", "mean")]

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







