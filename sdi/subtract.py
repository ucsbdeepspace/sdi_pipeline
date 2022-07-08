import click
import ois
#import zogy.py
from astropy.io.fits import CompImageHDU
from . import _cli as cli
from .combine import combine
from . import zogy_kyle_edit as zogy
import time

def subtract(hduls, name="ALGN", method: ("ois", "numpy", "zogy") ="ois"):
    """
    Returns differences of a set of images from a template image
    Arguments:
        hduls -- list of fits hdul's where the last image is the template 
        name -- name of the HDU to use image that the other images will be subtracted from
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    hduls = [h for h in hduls]
    outputs = []
    template = combine(hduls, name)["PRIMARY"].data # this is a temporary HDUL containing 1 PrimaryHDU with combined data ## Align is the reference image here??
    template_fits = combine(hduls, name)["PRIMARY"]
    print(" ")
    print("Beginning Subtraction")
    start = time.perf_counter()
    if method == "ois":
        print("Method = OIS")
        for i,hdu in enumerate(hduls):
            try:
                diff = ois.optimal_system(image=hdu[name].data, refimage=template, method='Bramich')[0] #Align is the input image here??
            except ValueError:
                diff = ois.optimal_system(image=hdu[name].data.byteswap().newbyteorder(), refimage=template.byteswap().newbyteorder(), method='Bramich')[0]
            hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i]['ALGN'].header, name = "SUB"))
            outputs.append(hdu)
    elif method == "numpy":
        print("Method = Numpy")
        for i,hdu in enumerate(hduls):
            diff = template - hdu[name].data
            hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i]['ALGN'].header, name = "SUB"))
            outputs.append(hdu)

    elif method == "zogy":
        print("Method = ZOGY")
        for i,hdu in enumerate(hduls):
            diff = zogy.optimal_subtraction(new_fits= hdu[name],      ref_fits=template_fits,
                        new_fits_mask=None, ref_fits_mask=None,
                        set_file='set_zogy', logfile=None,
                        redo_new=None, redo_ref=None,
                        verbose=None, nthreads=5, telescope='ML1',
                        keep_tmp=None)
            hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i]['ALGN'].header, name = "SUB"))
            outputs.append(hdu)
    else:
        raise ValueError(f"method {method} unknown!")
    stop = time.perf_counter()
    print("Subtraction Complete")
    print("Time Elapsed: {} sec".format(round(stop-start, 4)))
    return (hdul for hdul in outputs)

@cli.cli.command("subtract")
@click.option("-n", "--name", default="ALGN", help="The HDU to be aligned.")
@click.option("-m", "--method", default="ois", help="The subtraction method to use; ois or numpy (straight subtraction) or zogy.")
@cli.operator

## subtract function wrapper
def subtract_cmd(hduls, name="ALGN", method="ois"):
    """
    Returns differences of a set of images from a template image\n
    Arguments:\n
        hduls -- list of fits hdul's where the last image is the template\n
        name -- name of the HDU to use image that the other images will be subtracted from\n
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    return subtract(hduls, name, method)
