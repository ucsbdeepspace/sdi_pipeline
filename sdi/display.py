from . import _cli as cli
import click
from pyds9 import DS9
import numpy as np

def display(hduls, cats=None, color='green', size=40):
    """
    Opens a series of HDU Images on ds9, optionally with circled catalogs.
    :param hduls: list of HDUL to display
    :param cats: numpy recarray of catalog data to use. if it is a string or
        integer, it will be assumed that it refers to a specific BinTableHDU
        in the provided HDULists.
    :param color: color to circle things in
    :param size: size of circles
    """
    ds9 = DS9()
    ds9.set("scale mode zscale")
    # TODO Make this a real generator
    hduls = [h for h in hduls]
    try:
        cats = [hdul[cats].data for hdul in hduls]
    except (KeyError, AttributeError, IndexError):
        pass
    if cats is None:
        cats = (() for _ in hduls)

    for hdul, cat in zip(hduls, cats):
        ds9.set("frame new")
        ds9.set_pyfits(hdul)
        ds9.set("zoom to fit")
        for source in cat:
            ds9.set('regions command {{circle {} {} {} #color={}}}' \
                    .format(source['x'], source['y'], size, color))

    ds9.set("frame first")
    ds9.set("frame delete")
    
    return (hdul for hdul in hduls)

@cli.cli.command("display")
@click.option("-e", "--cat_ext", type=str, \
              help="EXTNAME for a TableHDU in each HDUL which "
              "contains a record array of sources to circle in ds9.")
@click.option("-c", "--color", default="green", help="A string color to draw"
              "circles in, e.g. 'green'")
@click.option("-s", "--size", default=40, help="The radius of the circles"
              "drawn", type=int)
@cli.operator
def display_cmd(hduls, cat_ext=None, color="green", size=40):
    """
    Opens multiple hduls on DS9 with the option of circling sources.
    from fits catalog tables.
    Continues to pass unmodified hduls down the chain
    """
    return display(hduls, cat_ext, color, size)
