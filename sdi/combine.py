"""
Combine merges a set of astronomical data from a template
History:
    Created/Extensively Refactored 2019-09-05
        Andrew Bluth <abluth@ucsb.edu>
"""
import click
from . import _cli as cli
 
import numpy as np
from astropy.io import fits

def combine(hduls, name="SCI"):
    """
    Combine takes a pixel-by-pixel median of a set of astronomical data to
    create a template image.
    Combine is a reduction. This means that the stream will be truncated and it
    will return just one image.

    \b
    :param hduls: list of fits hdul's
    :param name: the name of the HDU to sum among the HDULS
    :returns: a list with a single hdul representing the median image.
    """
    hduls_list = [hdul for hdul in hduls]
    try:
        data = [hdul[name].data for hdul in hduls_list]
    except KeyError as ke:
        hduls_list[0].info()
        raise KeyError(str(f"Name {name} not found in HDUList! Try running again with `combine -n [name]` from above")) from None
    comb = np.median(data, axis=0)
    hdu = fits.PrimaryHDU(comb)
    hduls_list += [fits.HDUList([hdu])]
    return fits.HDUList([hdu])

@cli.cli.command("combine")
@click.option("-n", "--name", default="SCI", help="The HDU to be aligned.")
@cli.operator
# TODO add option to pick out a specific table instead of just science

## combine function wrapper
def combine_cmd(hduls, name="SCI"):
    """
    Combine takes a pixel-by-pixel median of a set of astronomical data to
    create a template image.
    Combine is a reduction. This means that the stream will be truncated and it
    will return just one image.

    \b
    :param hduls: list of fits hdul's
    :param name: the name of the HDU to sum among the HDULS
    :returns: a list with a single hdul representing the median image.
    """
    return [combine(hduls, name),]
