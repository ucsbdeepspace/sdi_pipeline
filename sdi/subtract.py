import click
import ois
from astropy.io import fits                           #sfft specific
from astropy.io.fits import CompImageHDU
from . import _cli as cli
from .combine import combine
import os                                             #sfft specific
from sfft.EasyCrowdedPacket import Easy_CrowdedPacket #sfft specific
from sfft.CustomizedPacket import Customized_Packet   #sfft specific
from sfft.EasySparsePacket import Easy_SparsePacket   #sfft specific
import time
import numpy as np                                           #sfft specific
import sys

def subtract(hduls, name="ALGN", method: ("sfft", "ois", "numpy")="sfft"):
    """
    Returns differences of a set of images from a template image
    Arguments:
        hduls -- list of fits hdul's where the last image is the template 
        name -- name of the HDU to use image that the other images will be subtracted from
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    hduls = [h for h in hduls]
    outputs = []
    if method == "sfft":
        print(" ")
        print("Writing Temporary Fits Files")
        start = time.perf_counter()
        temp_image_fits_filenames = []
        venv_file_path = sys.prefix
        for i, hdu in enumerate(hduls): #Write the science hduls into temporary fits files
            temp_filename = venv_file_path + "/sdi_pipeline/sdi/temp_image{}.fits".format(i)
            temp_data = hdu[name].data
            primary = fits.PrimaryHDU(temp_data)
            primary.writeto(temp_filename, overwrite = True)
            temp_image_fits_filenames.append(temp_filename)
        temp_ref_path = venv_file_path + "/sdi_pipeline/temp_ref.fits"
        template_fits = combine(hduls, name) # Create the reference image
        template_fits.writeto(temp_ref_path, overwrite = True) # Write the reference image into temporary fits file

        stop = time.perf_counter()
        print("Writing Complete")
        print("Time Elapsed: {} sec".format(round(stop-start, 4)))
        
    # Check file size of temporary fits files on disc
        size = 0
        for i, fits_name in enumerate(temp_image_fits_filenames): 
            size += os.path.getsize(fits_name)
        size += os.path.getsize(temp_ref_path)
        print("Total size of temporary files is {} mb".format(size/1024**2))
        outputs = []

    #Subtract
        start = time.perf_counter()
        print("Method = sfft")  #TODO Incorperate input masks for when we take real data
        for i, fits_name in enumerate(temp_image_fits_filenames):       
            sol, diff = Customized_Packet.CP(FITS_REF = temp_ref_path, FITS_SCI = fits_name, 
                                    FITS_mREF = temp_ref_path, FITS_mSCI = fits_name,
                                    ForceConv = "REF", GKerHW = 4, BGPolyOrder = 2, KerPolyOrder = 2)    
            
            if np.isnan(np.sum(diff)) == True:
                raise ValueError("Residual contains NaN")
            else:
                hduls[i].insert(1,CompImageHDU(data = diff, header =  hduls[i][name].header, name = "SUB"))
            outputs.append(hduls[i])
        stop = time.perf_counter()
        print("Subtraction Complete")
        print("Time Elapsed: {} sec".format(round(stop-start, 4)))
        print("Removing Temporary Fits Files")    
        for i, fits_name in enumerate(temp_image_fits_filenames): #Remove Temporary Fits files from disc  
            if os.path.exists(fits_name):
                os.remove(fits_name)
            else:
                print("{} does not exist".format(fits_name))
        if os.path.exists(temp_ref_path):
                os.remove(temp_ref_path)
        else:
            print("temp_ref.fits does not exist")
        print("Removal Complete")
    elif method == "ois":
        print("Method = OIS")
        template = combine(hduls, name)
        for i,hdu in enumerate(hduls):
            try:
                diff = ois.optimal_system(image=hdu[name].data, refimage=template['PRIMARY'].data, method='Bramich')[0]
            except ValueError:
                diff = ois.optimal_system(image=hdu[name].data.byteswap().newbyteorder(), refimage=template.byteswap().newbyteorder(), method='Bramich')[0]
            hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i][name].header, name = "SUB"))
            outputs.append(hdu)

    elif method == "numpy":
        print("Method = Numpy")
        template = combine(hduls, name)
        for i,hdu in enumerate(hduls):
            diff = template["PRIMARY"].data - hdu[name].data
            hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i][name].header, name = "SUB"))
            outputs.append(hdu)
    else:
        raise ValueError(f"method {method} unknown!")
    return (hdul for hdul in outputs)

@cli.cli.command("subtract")
@click.option("-n", "--name", default="ALGN", help="The HDU to be aligned.")
@click.option("-m", "--method", default="sfft", help="The subtraction method to use; ois or numpy (straight subtraction) or sfft (GPU accelerated). sfft-sparse and sfft-croweded are in development")
@cli.operator

## subtract function wrapper
def subtract_cmd(hduls,name="ALGN", method="sfft"):
    """
    Returns differences of a set of images from a template image\n
    Arguments:\n
        hduls -- list of fits hdul's where the last image is the template\n
        name -- name of the HDU to use image that the other images will be subtracted from\n
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    return subtract(hduls, name, method)
