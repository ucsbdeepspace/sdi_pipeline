import os
import click
import ois
from astropy.io import fits
from . import _cli as cli
from .combine import combine

def subtract(hduls, name="SCI", method: ("ois", "numpy")="ois"):
    """
    Returns differences of a set of images from a template image
    Arguments:
        hduls -- list of fits hdul's where the last image is the template 
        name -- name of the HDU to use image that the other images will be subtracted from
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    hduls = [h for h in hduls]
    outputs = []
    template = combine(hduls, name)["PRIMARY"].data
    i = 0

    if method == "ois":
        for hdu in hduls:
            try:
                diff = ois.optimal_system(image=hdu[name].data, refimage=template, method='Bramich')[0]
            except:
                diff = ois.optimal_system(image=hdu[name].data.byteswap().newbyteorder(), refimage=template.byteswap().newbyteorder(), method='Bramich')[0]
            hdu.insert(0,fits.PrimaryHDU(diff))
            hdu.insert(1,fits.CompImageHDU(data = none, header =  hduls[i]['SCI'].header))
            outputs.append(hdu)
            i+=1

    elif method == "numpy":
        for hdu in hduls:
            diff = template - hdu[name].data
            hdu.insert(0,fits.PrimaryHDU(diff))
            hdu.insert(1,fits.CompImageHDU(data = none, header =  hduls[i]['SCI'].header))
            outputs.append(hdu)
            i+=1

    else:
        raise ValueError(f"method {method} unknown!")
    return (hdul for hdul in outputs)

@cli.cli.command("subtract")
@click.option("-n", "--name", default="SCI", help="The HDU to be aligned.")
@click.option("-m", "--method", default="ois", help="The subtraction method to use; ois or numpy (straight subtraction).")
@cli.operator

## subtract function wrapper
def subtract_cmd(hduls,name="SCI", method="ois"):
    """
    Returns differences of a set of images from a template image\n
    Arguments:\n
        hduls -- list of fits hdul's where the last image is the template\n
        name -- name of the HDU to use image that the other images will be subtracted from\n
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    return subtract(hduls, name, method)
