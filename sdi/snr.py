"""
snr calculates the signal-to-noise ratio of a fits file
Standalone version without click implementation
Code originally by Yael Brynjegard-Bialik 2020-2021
"""


import numpy as np
from astropy.io import fits
import os
from glob import glob
import sys as sys
from photutils import Background2D, DAOStarFinder, detect_threshold, detect_sources, source_properties
import matplotlib.pyplot as plt

def snr(hduls, name="SCI"):
	best = hduls[0]
	best_snr = 0
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

		source_max_values = SourceCatalog.max_value
		avg_source_max_values = np.mean(source_max_values)

		# calculate signal to noise ratio
		signal = avg_source_max_values
		noise = bkg_mean_rms
		SNR = (signal)/(noise)
		if SNR > best_snr:
			best = hdul
			best_snr = SNR
		hdul[name].header.append(('SNR',SNR,"signal to noise ratio" ))
	return best
