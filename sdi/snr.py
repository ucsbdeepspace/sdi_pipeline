"""
Standalone version without click implementation
Code originally by Yael Brynjegard-Bialik 2020-2021
Updated to return modified hduls instead of best SNR hdul by Tyler Pruitt 05-01-2021
"""


import numpy as np
from astropy.io import fits
import os
from glob import glob
import sys as sys
from photutils import Background2D, DAOStarFinder, detect_threshold, detect_sources, source_properties
import matplotlib.pyplot as plt

def snr(hduls, name="SCI"):
    """
    snr() calculates the signal-to-noise ratio (SNR) of a fits file and writes it to the header
    
    Parameters
    ----------
    hduls : list of hduls
        a connection of hduls without the SNR 
    
    name : string, optional
        name of the hdul. The default is "SCI".
    
    Returns
    -------
    hduls : list of hduls
        a modified connection of hduls with the SNR value written to the FITS header
    """
    
    for hdul in hduls:
   		data = hdul[name].data
   		shape = data.shape
   		
        # identify background rms
   		boxsize=(shape)
   		bkg = Background2D(data, boxsize)
   		bkg_mean_rms = np.mean(bkg.background_rms)
   
   		# subtract bkg from image
   		new_data = data - bkg.background
   
   		# set threshold and detect sources, threshold 5*std above background
   		threshold = detect_threshold(data=new_data, nsigma=5.0, background=0.0)
   		SegmentationImage = detect_sources(data=new_data, threshold=threshold, npixels=10)
   
   		SourceCatalog = source_properties(new_data, SegmentationImage)
   		columns = ['id', 'xcentroid', 'ycentroid', 'source_sum']
        
        #calculate mean max values of source
   		source_max_values = SourceCatalog.max_value
   		avg_source_max_values = np.mean(source_max_values)
   
   		# calculate signal to noise ratio
   		signal = avg_source_max_values
   		noise = bkg_mean_rms
   		SNR = (signal)/(noise)
        
        # write SNR value to the FITS header
   		hdul[name].header.append(('SNR', SNR, 'signal to noise ratio'))
    
    return (hdul for hdul in hduls)
