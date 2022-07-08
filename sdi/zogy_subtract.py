import click
import ois
import os
from astropy.io.fits import CompImageHDU
from . import _cli as cli
from .combine import combine
from . import zogy as zogy
import time

def zogy_subtract(hduls, name="ALGN"):
    """
    Returns differences of a set of images from a template image
    Arguments:
        hduls -- list of fits hdul's where the last image is the template 
        name -- name of the HDU to use image that the other images will be subtracted from
    """
    hduls = [h for h in hduls]
    # print("hduls info") #Remove before pushing
    # for i in hduls:
    #     print(i.info())
    print(" ")
    print("Writing Temporary Fits Files")
    start = time.perf_counter()
    temp_image_fits_filenames = []
    for i, hdu in enumerate(hduls):
        temp_filename = "temp_image{}.fits".format(i)
        hdu[name].writeto(temp_filename, overwrite = True)
        temp_image_fits_filenames.append(temp_filename)

    template_fits = combine(hduls, name)
    template_fits.writeto("temp_ref.fits", overwrite = True)

    stop = time.perf_counter()
    print("Writing Complete")
    print("Time Elapsed: {} sec".format(round(stop-start, 4)))
    size = 0
    for i, fits_name in enumerate(temp_image_fits_filenames):
        size += os.path.getsize(fits_name)
    size += os.path.getsize("temp_ref.fits")
    print("Total size of temporary files is {} mb".format(size/1024**2))
    outputs = []
    
    print(" ")
    print("Beginning Zogy Subtraction")
    start = time.perf_counter()
    print("Method = ZOGY")
    for i, fits_name in enumerate(temp_image_fits_filenames):
        diff = zogy.optimal_subtraction(new_fits= fits_name,      ref_fits="temp_ref_fits",
                new_fits_mask=None, ref_fits_mask=None,
                set_file='set_zogy', logfile=None,
                redo_new=None, redo_ref=None,
                verbose=None, nthreads=5, telescope='ML1',
                keep_tmp=None)
    hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i]['ALGN'].header, name = "SUB"))
    outputs.append(hdu)
    stop = time.perf_counter()
    print("Subtraction Complete")
    print("Time Elapsed: {} sec".format(round(stop-start, 4)))
    print("Removing Temporary Fits Files")
    for i, fits_name in temp_image_fits_filenames:
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

@cli.cli.command("zogy_subtract")
@click.option("-n", "--name", default="ALGN", help="The HDU to be aligned.")
@cli.operator

## subtract function wrapper
def zogy_subtract_cmd(hduls, name="ALGN"):
    """
    Returns differences of a set of images from a template image\n
    Arguments:\n
        hduls -- list of fits hdul's where the last image is the template\n
        name -- name of the HDU to use image that the other images will be subtracted from\n
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    return zogy_subtract(hduls, name)
