# Preface: This "runscript" is branched from "runscript_serial_membound.py". 
# Whereas "runscript_serial_membound.py" tried to solve the memory issue
# that arises during spherematching in Step 4 by dividing the computation
# in chunks and accumulate the results in a list, the present script takes
# in two command line arguments that specify which "chunk" of the HEALPix
# pixels to process. The two arguments specify the start and end index and
# used for slicing. See the actual code for the description of the workflow.
# For higher level explanation of steps, see below.
# 
# Use "config_run_per_chunk.py" file as the configuration file.
# 
# Given the config file provided by the user, the present script produces the
# output array of interest. Here is a high-level overview of the program.
# 
# - 1. Check config.py: The runscript checks whether the config file is properly 
#   written and produces any warning signs that the user should be aware of. 
#   It aborts the program if the config file is not proper. The file should 
#	be self-explanatory.
# 
# - 2. Load ccd data: Load the ccd summary file located in the dirctory
#	specified in the config file. Mask out rows according to the conditions.
#	This part is possibly subject to review. Extract the columns the user is 
# 	interested and preprocess them according to the rules specified by the 
#	functions.py. Return trimmed ccd data and in particular its ra/dec.
# 
# - 3. Create HEALPix grid: Using Nside specified by the user generate HEALPix
# 	grid with nest=NESTED, lonlat=True. Default Nside = 2^11 is recommended
#	for speed and accuracy of statistics computed for each pixel.
#
# - 4. Spherematch HEALPix centers to ccd centers within 0.336/2 degrees, half the
#	size of the diagonal of ccd frame. Use astropy.coordinates.search_around_sky.
#	The output is idx_pix, idx_ccd that gives indices of all matches.
# 	`
# - 5. Trim the list of matches: For each HEALPix pixel (referred to as pixel 
#	from here on), find a set of ccds that it actually belongs to. This can 
#	be done efficiently in cartesian represntation of ra/dec's as described 
#	in the code below. The cartesian represntation is pre-computed. 
#
# - 6. Compute the stats: The statistics specified by the user (e.g. min, max, mean)
#	are computed via scipy.stats.binned_statistic function.
#
# - 7. Report times for various steps. 
#
# - Plots: Currently not supported.
#
# TODO: When the program runs, first check whether the inputs are reasonable.
# TODO: What to do with galdepth or psfdepth with zero? Or if the values don't make sense?
# TODO: If the user specifies a silly template operation, then print WARNING and sometimes abort.



################################################################################
# Required modules/functions
import healpy as hp
import numpy as np
import matplotlib.pylab as plt
import astropy.io.fits as fits 
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
from astropy import units as u
import scipy.stats as stats
import time
from multiprocessing import Pool
# module of ccdtosky
from functions import *
# read in configuration file
from config_run_per_chunk import *
# For reading in command line args
import sys



################################################################################
# - 1. Check config.py: The runscript checks whether the config file is properly 
#   written and produces any warning signs that the user should be aware of. 
#   It aborts the program if the config file is not proper. The file should 
#	be self-explanatory.
print("1. Check config.py")
start = time.time()
dt1 = time.time()-start
# Take the command line inputs for start and end indices for the HEALPix pixel
# chunk to compute.
chunk_start = int(sys.argv[1])
chunk_end = int(sys.argv[2])
num_pix_tot = hp.nside2npix(Nside)
# If the chunk_end exceeds the bound, then set it the last possible value.
if chunk_end > num_pix_tot:
    chunk_end = num_pix_tot
# If start is greater than end, then abort.
if chunk_end <= chunk_start:
    raise Exception("Error: Chunk start index is greater than end index.")
print("Finished. Time elapsed: %.3E sec\n"% (dt1))
print("\n")



################################################################################
# - 2. Load ccd data: Load the ccd summary file located in the dirctory
#	specified in the config file. Mask out rows according to the conditions.
#	This part is possibly subject to review. Extract the columns the user is 
# 	interested and preprocess them according to the rules specified in the 
#	functions.py. Return trimmed ccd data and in particular its ra/dec.
print("2. Load ccd data")
start = time.time()
data_ccd= load_ccd_data(faddress)
num_ccds_total = data_ccd.size
print("ccd file loaded is loaded from %s." % faddress)
print("Total # ccds in the summary file: {:>,d}".format(num_ccds_total))

# In the following, select rows corresponding to ccd frames that were used in data reduction.
# There are three possible conditions that would exclude a particular ccd:
# 1) "photometric"
# 2) "blacklist"
# 3) "goodregion"
# Of the three, I found that only "photometric" and "good region" condition imposes any restriction.
# To make the matter simple, for now, I throw away any frames that do not pass all three conditions.

# Compute the number of frames with various conditions. (Code intentionally verbose for readability.)
iphoto = data_ccd["photometric"]==True
iblacklist_ok = data_ccd["blacklist_ok"]==True
igood_region = np.sum(data_ccd["good_region"],axis=-1)==-4

num_ccd_photo = np.sum(iphoto)
num_ccd_blacklist_ok = np.sum(iblacklist_ok)
num_ccd_good_region = np.sum(igood_region)

# Report the computed numbers above.
print("Masking ccd frames with poor conditions.")
print("# ccds with \"photometric\" bad: {:>,d} ({:2.2f} pcnt)".format(num_ccds_total- \
            num_ccd_photo,(num_ccds_total-num_ccd_photo)/num_ccds_total * 100))
print("# ccds blacklisted: {:>,d} ({:2.2f} pcnt)".format(num_ccds_total- \
            num_ccd_blacklist_ok,(num_ccds_total-num_ccd_blacklist_ok)/num_ccds_total * 100))
print("# ccds only partially good: {:>,d} ({:2.2f} pcnt)".format(num_ccds_total- \
            num_ccd_good_region,(num_ccds_total-num_ccd_good_region)/num_ccds_total * 100))

# Boolean vector for the intersection of the three conditions.
ibool = np.logical_and.reduce((iphoto, iblacklist_ok, igood_region))
num_ccd_used = np.sum(ibool)
print("# ccds USED after masking: {:>,d} ({:2.2f} pcnt)".format(num_ccd_used,(num_ccd_used)/num_ccds_total * 100))

# Overwritting data_ccd variable so as to use only the unmasked data.
data_ccd = data_ccd[ibool]

# The centers of ccd
# print("\n")
ra_ccd, dec_ccd = load_radec_center(data_ccd)
print("ccd ra/dec centers loaded.")
dt2 = time.time()-start
print("Finished. Time elapsed: %.3E sec\n"% (dt2))
print("\n")



################################################################################
# - 3. Create HEALPix grid: Using Nside specified by the user generate HEALPix
# 	grid with nest=NESTED, lonlat=True. Default Nside = 2^11 is recommended
#	for speed and accuracy of statistics computed for each pixel.
print("3. Create HEALPix grid")
start = time.time()
print("Nside chosen: %d"% Nside)
print("Total # HEALPix pixels: %d"%hp.nside2npix(Nside))
num_pix = -(chunk_start-chunk_end) 
print("# pixels computed: %d"%num_pix)
print("Computing chunk [%d, %d]" % (chunk_start, chunk_end))
ra_pix, dec_pix = np.array(hp.pix2ang(Nside,range(chunk_start, chunk_end), nest=NESTED, lonlat=True)) # lonlat=True, must.
print("HEALPix pixels ra/dec computed.")

### For each healpix pixel position, finding ccd centers.
c_pix = SkyCoord(ra=ra_pix*u.degree, dec=dec_pix*u.degree)
c_ccd = SkyCoord(ra=ra_ccd*u.degree, dec=dec_ccd*u.degree)
print("Astropy SkyCoord objs for ccd and pix ra/dec created.")
dt3 = time.time()-start
print("Finished. Time elapsed: %.3E sec\n"% (dt3))
print("\n")



################################################################################
# - 4. Spherematch HEALPix centers to ccd centers within 0.336/2 degrees, half the
#	size of the diagonal of ccd frame. Use astropy.coordinates.search_around_sky.
#	The output is idx_pix, idx_ccd that gives indices of all matches.
print("4. Spherematch HEALPix centers to ccd centers within 0.336/2 degrees")
print("Computing ccd to pix mapping based on their ra/dec's. (kD-tree)")
start = time.time()

# We make c_pix the first argument because we want the mapping to be sorted by it.
idx_pix, idx_ccd, _, _ = search_around_sky(c_pix, c_ccd, seplimit=sepdeg*u.degree, storekdtree='kdtree_sky')

dt4 = time.time()-start 
print("Finished. Time took: %.3E sec\n" %(dt4))

# If the number of matches is zeros, then exit the program.
if idx_pix.size == 0:
    print("Number of matches (length of idx_pix): {:>,d}".format(idx_pix.size))
    print("Total # HEALpix pixels: %d" % num_pix_tot)
    print("# pixels computed: %d"%num_pix)    
    print("Computing chunk [%d, %d]" % (chunk_start, chunk_end))
    print("{:<40s}: {:<1.3E} sec".format("Check config file", dt1))
    print("{:<40s}: {:<1.3E} sec".format("Load CCD data", dt2))
    print("{:<40s}: {:<1.3E} sec".format("Create HEALPix grid", dt3))
    print("{:<40s}: {:<1.3E} sec".format("Spherematch pixels to ccd centers", dt4))
    print("# ccds USED after masking: {:>,d} ({:2.2f} pcnt)".format(num_ccd_used,(num_ccd_used)/num_ccds_total * 100))
    print("Skip the rest of the program and exit.")
    sys.exit()

print("In the remainder, work with only pixels that have matching ccd frames.")
idx_pix_uniq = np.unique(idx_pix) # Get the uniq set of pixels matchted.
num_pix_uniq = idx_pix_uniq.size
print("# pix that overlaps with imaged regions: %d (%2.2f pcnt)" %(num_pix_uniq, num_pix_uniq/num_pix*100))
print("Number of matches (length of idx_pix): {:>,d}".format(idx_pix.size))
print("\n")



################################################################################	
# - 5. Trim the list of matches: For each HEALPix pixel (referred to as pixel 
#	from here on), find a set of ccds that it actually belongs to. This can 
#	be done efficiently in cartesian represntation of ra/dec's as described 
#	in the code below. The cartesian represntation is pre-computed. 
print("5. Trim the list of matches")
print("5a: Loading ra/dec bounding corners and covert them (and pixel ra/dec) to xyz on unit sphere.")
start = time.time()

# Load ra/dec of ccd bounding corners.
ra0, ra1, ra2, ra3, dec0, dec1, dec2, dec3 = load_radec_corners(data_ccd)

# Convert ra/dec to xyz vectors
xyz0_ccd = radec2xyz(ra0,dec0)
xyz1_ccd = radec2xyz(ra1,dec1)
xyz2_ccd = radec2xyz(ra2,dec2)
xyz3_ccd = radec2xyz(ra3,dec3)

# Compute the normal vectors. The idea here is that given four normal vectors that 
# are the results of cross-producting (xyz1,xyz0), (0,3), (3,2), (2,1) in that order, 
# we can simply dot prodcut with normal vectors represented by healpix center.
# The subtraction is done for numerically stability. Members of n0_ccd-n_pix 
# are quite similar.
n0_ccd  = np.cross(xyz1_ccd-xyz0_ccd, xyz0_ccd)
n1_ccd  = np.cross(xyz0_ccd-xyz3_ccd, xyz3_ccd)
n2_ccd  = np.cross(xyz3_ccd-xyz2_ccd, xyz2_ccd)
n3_ccd  = np.cross(xyz2_ccd-xyz1_ccd, xyz1_ccd)

# Create x,y,z vector version of healpix centers.
xyz_pix = np.asarray(c_pix.cartesian.xyz).T
dt5a = time.time()-start
print("Finished. Time elapsed: %.3E sec\n"% (dt5a))

# Vectorized version. 
print("5b. Trimming ccd to pix mapping to only pairs where pix center resides inside the matching ccd.")
start = time.time()

# Get xyz vector of healpix center
n_pix = xyz_pix[idx_pix]

# Test whether the pixel center is within each of the matched CCD
ibool = np.logical_and.reduce((vectorized_dot(n0_ccd[idx_ccd],n_pix)>0, 
                               vectorized_dot(n1_ccd[idx_ccd],n_pix)>0, 
                               vectorized_dot(n2_ccd[idx_ccd],n_pix)>0, 
                               vectorized_dot(n3_ccd[idx_ccd],n_pix)>0))

# Trimming the mapping.
idx_pix_inside = idx_pix[ibool]
idx_ccd_inside = idx_ccd[ibool]
dt5b = time.time()-start
print("Finished. Time elapsed: %.3E sec\n"% (dt5b))

print("From here on work with pixels that were found to be in at least one ccd.")
idx_pix_inside_uniq = np.unique(idx_pix_inside)
num_pix_inside_uniq = idx_pix_inside_uniq.size
print("Pix # Beginning, # Spherematched, # Inside: %d, %d, %d "%(num_pix,num_pix_uniq,num_pix_inside_uniq))
print("\n")



################################################################################
# - 6. Compute the stats: See above.
print("6. Compute the stats")
print("Allocate memory for the output array.")
# Create the output recarray as a placeholder.
rec_dtype = gen_rec_dtype(templates)
num_col = 3*len(templates)+3
output_arr = np.recarray((num_pix_inside_uniq,),dtype=rec_dtype)

# Overwrite for the HEALPix portion.
output_arr["hpix_idx"] = idx_pix_inside_uniq
output_arr["hpix_ra"],output_arr["hpix_dec"] = hp.pix2ang(Nside,idx_pix_inside_uniq,nest=True,lonlat=True)

# Filter types
filter_types = ["g","r","z"]

# Before proceeding, if any of the templates specified "weight" = "galdepth_ivar", then compute the quanity.
for e in templates:
    if e[1]=="galdepth_ivar":
        print("Compute galdepth_ivar.")
        galdepth_ivar = fluxlim2ivar(mag2flux(data_ccd["galdepth"]))
        break

# Import ccd_filter
ccd_filter = data_ccd["filter"]

# Look-up dictionary for efficient assigning of numbers after computation.
eb_dict = template_filter_dict(templates, filter_types)


# Version: Using scipy binned_statistic
print("Start computing statistics.")
start = time.time()

# For each filter
for b in filter_types:
    print("Computation for %s-band quantities started."%b)	
    # Create a band mask
    i_b = data_ccd["filter"][idx_ccd_inside]==b

    # Sum of galdepth_ivar per bin
    # For galdepth_ivar average:
    hist_denom_galdepth_ivar, _, _ = stats.binned_statistic(idx_pix_inside[i_b], galdepth_ivar[idx_ccd_inside[i_b]], statistic = "sum", bins=np.arange(-0.5, num_pix+1.5, 1))

    # For each quantity requsted
    for e in templates:
        start_e = time.time()
        if (e[1] == "none") or (e[2] in ["min", "max", "median"]): 
            # If the weight scheme is none or any of the above functions were chosen.
            weight = False
        else:
            weight = True
        
        # If the quantity requested is Nexp
        if (e[0] == "Nexp"): # Nexp is not part of ccd file summary so treated like a special case
            hist_num, _, _= stats.binned_statistic(idx_pix_inside[i_b], np.ones_like([idx_ccd_inside[i_b]]), statistic = "sum", bins=np.arange(-0.5, num_pix+1.5, 1))            
            output_arr[eb_dict[(e,b)]] = hist_num[0][idx_pix_inside_uniq]            
        else:
            # If the operation asked for is mean
            if e[2] == "mean": 
                if weight: # weight is true, apply the weights when computing the average.
                    hist_num, _, _= stats.binned_statistic(idx_pix_inside[i_b], data_ccd[e[0]][idx_ccd_inside[i_b]]*galdepth_ivar[idx_ccd_inside[i_b]], statistic = "sum", bins=np.arange(-0.5, num_pix+1.5, 1))
                    output_arr[eb_dict[(e,b)]] = (hist_num[idx_pix_inside_uniq]/hist_denom_galdepth_ivar[idx_pix_inside_uniq])  
                else: # weight is FALSE, then just use "mean" option. 
                    hist_num, _, _= stats.binned_statistic(idx_pix_inside[i_b], data_ccd[e[0]][idx_ccd_inside[i_b]], statistic = "mean", bins=np.arange(-0.5, num_pix+1.5, 1))
                    output_arr[eb_dict[(e,b)]] = hist_num[idx_pix_inside_uniq]
            # For all other operations.
            else:
                hist_num, _, _= stats.binned_statistic(idx_pix_inside[i_b], data_ccd[e[0]][idx_ccd_inside[i_b]], statistic = e[2], bins=np.arange(-0.5, num_pix+1.5, 1))
                output_arr[eb_dict[(e,b)]] = hist_num[idx_pix_inside_uniq]
        print(("Quantity: %s, Time taken: {:<1.3E} sec"%eb_dict[(e,b)]).format(time.time()-start_e))
    print("Computation for %s-band quantities ended.\n"%b)    
                
dt6 = time.time()-start
print("Finished. Time elapsed: %.3E sec"% (dt6))
# At the moment, saved in numpy binary file.
np.save("".join([out_directory, "output_arr_chunk%dthru%d"%(chunk_start, chunk_end)]), output_arr, allow_pickle=False)
print("output_arr is saved in %s." % out_directory)
print("\n")



################################################################################
# - 7. Report times of various steps.
print("Total # HEALpix pixels: %d" % num_pix_tot)
print("# pixels computed: %d"%num_pix)
print("Computing chunk [%d, %d]" % (chunk_start, chunk_end))
print("{:<40s}: {:<1.3E} sec".format("Check config file", dt1))
print("{:<40s}: {:<1.3E} sec".format("Load CCD data", dt2))
print("{:<40s}: {:<1.3E} sec".format("Create HEALPix grid", dt3))
print("{:<40s}: {:<1.3E} sec".format("Spherematch pixels to ccd centers", dt4))
print("{:<40s}: {:<1.3E} sec".format("Pre-compute xyz coordinates", dt5a))
print("{:<40s}: {:<1.3E} sec".format("Trim the spherematch list", dt5b))
print("{:<40s}: {:<1.3E} sec".format("Compute the statistics", dt6))
print("# ccds USED after masking: {:>,d} ({:2.2f} pcnt)".format(num_ccd_used,(num_ccd_used)/num_ccds_total * 100))
print("Pix # Beginning, # Spherematched, # Inside: %d, %d, %d "%(num_pix,num_pix_uniq,num_pix_inside_uniq))
print("Number of matches (length of idx_pix): {:>,d}".format(idx_pix.size))




################################################################################
# - Plots: Currently not supported.