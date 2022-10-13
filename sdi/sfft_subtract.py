import click
import ois
import os
from astropy.io import fits
from astropy.io.fits import CompImageHDU
from . import _cli as cli
from .combine import combine
from sfft.EasyCrowdedPacket import Easy_CrowdedPacket
from sfft.CustomizedPacket import Customized_Packet
from sfft.EasySparsePacket import Easy_SparsePacket
import time
import numpy as np

def sfft_subtract(hduls, name="ALGN", method: ("Crowded, Sparse, None") = "None"):
    """
    Returns differences of a set of images from a template image
    Arguments:
        hduls -- list of fits hdul's where the last image is the template 
        name -- name of the HDU to use image that the other images will be subtracted from
        method -- whether to run image-masking preprocessing for sparse fields or crowded fields or no image-masking preprocessing
    """
    hduls = [h for h in hduls]

# Write temporary fits files to feed into SFFT subtract
# TODO: Remove if SFFT subtract updates to accept data stacks instead of full fits files
    print(" ")
    print("Writing Temporary Fits Files")
    start = time.perf_counter()
    temp_image_fits_filenames = []
    for i, hdu in enumerate(hduls): #Write the science hduls into temporary fits files
        temp_filename = "temp_image{}.fits".format(i)
        hdu[name].writeto(temp_filename, overwrite = True)
        temp_image_fits_filenames.append(temp_filename)

    template_fits = combine(hduls, name) # Create the reference image
    template_fits.writeto("temp_ref.fits", overwrite = True) # Write the reference image into temporary fits file

    stop = time.perf_counter()
    print("Writing Complete")
    print("Time Elapsed: {} sec".format(round(stop-start, 4)))
    
# Check file size of temporary fits files on disc
    size = 0
    for i, fits_name in enumerate(temp_image_fits_filenames): 
        size += os.path.getsize(fits_name)
    size += os.path.getsize("temp_ref.fits")
    print("Total size of temporary files is {} mb".format(size/1024**2))
    outputs = []

#Subtract
    print(" ")
    print("Beginning SFFT Subtraction")
    start = time.perf_counter()
    print("Method = SFFT")
    if method == "None":
        #TODO Incorperate input masks for when we take real data
        for i, fits_name in enumerate(temp_image_fits_filenames):       
            sol, diff = Customized_Packet.CP(FITS_REF = "temp_ref.fits", FITS_SCI = fits_name, 
                                    FITS_mREF = 'temp_ref.fits', FITS_mSCI = fits_name,
                                    ForceConv = "REF", GKerHW = 4, BGPolyOrder = 2, KerPolyOrder = 2)    
            hdul = fits.open(fits_name)
            if np.isnan(np.sum(diff)) == True:
                raise ValueError("Residual contains NaN")
            else:
                hdul.insert(1,CompImageHDU(data = diff, header =  hdul['ALGN'].header, name = "SUB"))
            outputs.append(hdul)
    #TODO Make Crowded Work for our Files
    elif method == 'Crowded':
        for i, fits_name in enumerate(temp_image_fits_filenames):
            prep, sol, diff = Easy_CrowdedPacket.ECP(FITS_REF = "temp_ref.fits", FITS_SCI = fits_name)
            hdul = fits.open(fits_name)
            hdul.insert(1,CompImageHDU(data = diff, header =  hdul['ALGN'].header, name = "SUB"))
            outputs.append(hdul)

    #TODO Make Sparse Work for our Files
    elif method == "Sparse":
        for i, fits_name in enumerate(temp_image_fits_filenames):
            prep, sol, diff = Easy_SparsePacket.ESP(FITS_REF = "temp_ref.fits", FITS_SCI = fits_name)
            hdul = fits.open(fits_name)
            hdul.insert(1,CompImageHDU(data = diff, header =  hdul['ALGN'].header, name = "SUB"))
            outputs.append(hdul)    
    else:
        raise TypeError("Method must either be 'Crowded', 'Sparse', or 'None'")    
    stop = time.perf_counter()
    print("Subtraction Complete")
    print("Time Elapsed: {} sec".format(round(stop-start, 4)))
    
#Remove Temporary Fits files from disc    
    print("Removing Temporary Fits Files")
    for i, fits_name in enumerate(temp_image_fits_filenames):
        if os.path.exists(fits_name):
            os.remove(fits_name)
        else:
            print("{} does not exist".format(fits_name))
    if os.path.exists("temp_ref.fits"):
            os.remove("temp_ref.fits")
    else:
        print("temp_ref.fits does not exist")
    print("Removal Complete")
    return (hdul for hdul in outputs)

@cli.cli.command("sfft_subtract")
@click.option("-n", "--name", default="ALGN", help="The HDU to be aligned.")
@click.option("-m", "--method", default = "None", help = "The type of image preprocessing to be used")
@cli.operator

## subtract function wrapper
def sfft_subtract_cmd(hduls, name="ALGN", method: ("Crowded", "Sparse", "None")= "None"):
    """
    Returns differences of a set of images from a template image\n
    Arguments:\n
        hduls -- list of fits hdul's where the last image is the template\n
        name -- name of the HDU to use image that the other images will be subtracted from\n
        method -- whether to run image-masking preprocessing for sparse fields or crowded fields
    """
    return sfft_subtract(hduls, name, method)
