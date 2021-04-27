import os
import click
import ois
from astropy.io import fits
import sdi
from .combine import combine

@sdi.cli.command("subtract")
@click.option("-n", "--name", default="SCI", help="The HDU to be aligned.")
@sdi.operator
def subtract(hduls, name="SCI"):
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
     
    for hdu in hduls:
        diff = ois.optimal_system(image=hdu[name].data, refimage=template, method='Bramich')[0]
        output.append(diff)     

    for item in output:
        # FIXME this is ragingly wrong, multiple items should be associated
        hdu = fits.PrimaryHDU(item)
        outputs.append(fits.HDUList([hdu])) 
    return (hdul for hdul in outputs)
