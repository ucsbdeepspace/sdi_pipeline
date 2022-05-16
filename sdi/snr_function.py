"""
snr -- this module calculates the signal-to-noise ratio of a fits file
Code originally by Yael Brynjegard-Bialik 2020-2021
"""


import numpy as np
from photutils import Background2D, detect_threshold, detect_sources, source_properties
import click
from . import _cli as cli
import time
def snr(hduls, name="SCI"):
    """
    calculates the signal-to-noise ratio of a fits file
    """
    time_start = time.time()
    times = []
    for hdul in hduls:
        loop_start = time.time()
        data = hdul[name].data

        # identify background rms
        boxsize = (data.shape)
        bkg = Background2D(data, boxsize)
        bkg_mean_rms = np.mean(bkg.background_rms)

        # subtract bkg from image
        new_data = data - bkg.background

        # set threshold and detect sources, threshold 5*std above background
        threshold = detect_threshold(data=new_data, nsigma=5.0, background=0.0)
        segmentation_image = detect_sources(data=new_data, threshold=threshold, npixels=10)

        source_catalog = source_properties(new_data, segmentation_image)
        # columns = ['id', 'xcentroid', 'ycentroid', 'source_sum']

        source_max_values = source_catalog.max_value
        avg_source_max_values = np.mean(source_max_values)

        # calculate signal to noise ratio
        signal = avg_source_max_values
        noise = bkg_mean_rms
        sig_to_noise = (signal)/(noise)
        hdul["CAT"].header.append(('SNR', sig_to_noise, "signal to noise ratio"))
        times.append(time.time() - loop_start)
    print(f"AVG time to calc SNR : {np.mean(times)}")
    print(f"Time to run snr_function.py : {time.time() - time_start}")
    return (hdul for hdul in hduls)

@cli.cli.command("snr")
@click.option("-n", "--name", default="SCI", help="The HDU to calculate for")
@cli.operator

def snr_cmd(hduls, name="SCI"):
    """
    snr function wrapper
    """
    return snr(hduls, name)
