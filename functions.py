# Read variables from config. This shouldn't be necessary so I am commenting 
# it out.
# from config import *

# pakcages and modules
import healpy as hp
import numpy as np
import matplotlib.pylab as plt
import astropy.io.fits as fits 
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
from astropy import units as u
import time



def load_ccd_data(faddress="ccds-annotated-decals.fits"):
    """
    Load ccd data file given the file address.
    """
    return fits.open(faddress)[1].data


def load_radec_center(data_ccd):
    """
    Given the ccd data rec array, return ccd (ra,dec) cetners
    """
    return data_ccd["ra_center"],data_ccd["dec_center"]


def load_radec_corners(data_ccd):
    """
    Given the ccd data rec array, return ccd (ra,dec) cetners
    """
    return data_ccd["ra0"],data_ccd["ra1"],data_ccd["ra2"],data_ccd["ra3"],data_ccd["dec0"],data_ccd["dec1"],data_ccd["dec2"],data_ccd["dec3"]


def radec2xyz(ra,dec):
    """
    Covert ra/dec to xyz on a unit sphere.
    """
    return np.asarray(SkyCoord(ra=ra*u.degree, dec=dec*u.degree).cartesian.xyz).T


def mag2flux(mag):
    """
    AB magnitude to nanomaggie conversion.
    """
    return 10**((22.5-mag)/2.5)


def fluxlim2ivar(fluxlim, sig=5):
    """
    SNR_lim = (f_lim/df) = sig --> ivar = (1/df)^2 = (sig/f_lim)^2
    """
    return ((sig+1e-12)/(fluxlim).astype(np.float))**2


def gen_rec_dtype(templates):
    """
    Given the templates, generate dtype array for the output recarry.
    """
    filter_types = ["g","r","z"]
    rec_dtype = [("hpix_idx",int), ("hpix_ra",float), ("hpix_dec",float)]
    for e in templates:
        for b in filter_types:
            if callable(e[2]): # If the user passed the function
                if e[1]=="none":
                    rec_dtype.append(("_".join([b,e[0],"userfunc"]),float))
                else:
                    rec_dtype.append(("_".join([b,e[0],"w","userfunc"]),float))                    
            else:
                if e[1]=="none":
                    rec_dtype.append(("_".join([b,e[0],e[2]]),float))
                else:
                    rec_dtype.append(("_".join([b,e[0],"w",e[2]]),float))                    
    return rec_dtype


def template_filter_dict(templates, filter_types):
    """
    Returns a dictionary whose key is a tuple of (template element, filter)
    and the corresponding value is column name of the output array.
    """
    eb_dict = {}
    for e in templates:
        for b in filter_types:
            if callable(e[2]): # If the user passed the function
                ptr = "_".join([b,e[0],"userfunc"])
            else:
                ptr = "_".join([b,e[0],e[2]])
            eb_dict[(e,b)] = ptr
    return eb_dict

def compute_stat(template_element, data_ccd, galdepth_ivar, idx_ccd_matched):
    """
    Compute the statistic of the quantity specified in template_element, given ccd data and weights.
    Recall that the user could specify any of the following functions.
        - "mean": If "weight" = "galdepth_ivar", then weighted average is computed.
        - "min"
        - "max"
        - "median"
        - "sum"
        - "std": numpy.std.
    """
    # Unpack template element
    ccd_quantity, weight, func = template_element
    
    # Compute the quantity
    if ccd_quantity == "Nexp":
        qs = np.ones(idx_ccd_matched.size)
    else:
        qs = data_ccd[ccd_quantity][idx_ccd_matched] # Quantities
    
    # Compute the weight
    if weight == "none":
        ws = np.ones_like(qs)
    else:
        ws = galdepth_ivar[idx_ccd_matched] # Weights
    
    # Assign the proper function.
    if callable(func):
        ans = func(qs, ws)
    else:
        if func == "mean": ans = np.average(qs, weights=ws)
        if func == "min": ans = np.min(qs)
        if func == "max": ans = np.max(qs)
        if func == "median": ans = np.median(qs)
        if func == "sum": ans = np.sum(qs)
        if func == "std": ans = np.std(qs)
            
    return ans
    

def vectorized_dot(a,b):
    """
    Given two matrices of Nx3 where N represents number of vectors
    and 3 represents xyz components, return a vector of N giving dot 
    product corresponding to each row.
    """
    return np.sum(a*b,axis=1)