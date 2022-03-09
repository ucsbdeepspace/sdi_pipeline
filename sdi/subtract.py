import click
import ois
from astropy.io.fits import CompImageHDU
from . import _cli as cli
from .combine import combine

def subtract(hduls, name="ALGN", method: ("ois", "numpy")="ois"):
    """
    Returns differences of a set of images from a template image
    Arguments:
        hduls -- list of fits hdul's where the last image is the template 
        name -- name of the HDU to use image that the other images will be subtracted from
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    hduls = [h for h in hduls]
    outputs = []
    template = combine(hduls, name)["PRIMARY"].data # this is a temporary HDUL containing 1 PrimaryHDU with combined data
    if method == "ois":
        for i,hdu in enumerate(hduls):
            try:
                diff = ois.optimal_system(image=hdu[name].data, refimage=template, method='Bramich')[0]
            except ValueError:
                diff = ois.optimal_system(image=hdu[name].data.byteswap().newbyteorder(), refimage=template.byteswap().newbyteorder(), method='Bramich')[0]
            hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i]['ALGN'].header, name = "SUB"))
            outputs.append(hdu)

    elif method == "numpy":
        for i,hdu in enumerate(hduls):
            diff = template - hdu[name].data
            hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i]['ALGN'].header, name = "SUB"))
            outputs.append(hdu)
    else:
        raise ValueError(f"method {method} unknown!")
    return (hdul for hdul in outputs)

@cli.cli.command("subtract")
@click.option("-n", "--name", default="ALGN", help="The HDU to be aligned.")
@click.option("-m", "--method", default="ois", help="The subtraction method to use; ois or numpy (straight subtraction).")
@cli.operator

## subtract function wrapper
def subtract_cmd(hduls,name="ALGN", method="ois"):
    """
    Returns differences of a set of images from a template image\n
    Arguments:\n
        hduls -- list of fits hdul's where the last image is the template\n
        name -- name of the HDU to use image that the other images will be subtracted from\n
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    return subtract(hduls, name, method)
