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
        name -- name of the HDU to use
        image that the other images will be subtracted from
    """
    hduls = [h for h in hduls]
    output = []
    outputs = []
    template = combine(hduls, name)["PRIMARY"].data

    if method == "ois":
        for hdu in hduls:
            diff = ois.optimal_system(image=hdu[name].data.byteswap().newbyteorder(), refimage=template.byteswap().newbyteorder(), method='Bramich')[0]
            output.append(diff)     
    elif method == "numpy":
        for hdu in hduls:
            diff = template - hdu[name].data
            output.append(diff)     
    else:
        raise ValueError(f"method {method} unknown!")

    for item in output:
        newhdu = fits.HDUList([fits.PrimaryHDU(item)] + hdu[1:])
        outputs.append(newhdu)
    return (hdul for hdul in outputs)

@cli.cli.command("subtract")
@click.option("-n", "--name", default="SCI", help="The HDU to be aligned.")
@cli.operator

## subtract function wrapper
def subtract_cmd(hduls,name="SCI"):
    """
    Returns differences of a set of images from a template image
    Arguments:
        hduls -- list of fits hdul's where the last image is the template
        name -- name of the HDU to use
        image that the other images will be subtracted from
    """
    return subtract(hduls, name)
