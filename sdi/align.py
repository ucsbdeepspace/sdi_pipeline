"""
align -- this module aligns sets of astronomical data
HISTORY
    Created/Substantially refactored on 2019-09-01
        Varun Iyer <varun_iyer@ucsb.edu>
"""
# general imports

import click
import numpy as np
from astropy.io.fits import CompImageHDU
#from sources import Source
import astroalign
from _scripts import snr
from . import _cli as cli

def align(hduls, name="SCI", reference=None):
    """
    Aligns the source astronomical image(s) to the reference astronomical image
    \b
    :param hduls: list of fitsfiles
    :return: list of fistfiles with <name> HDU aligned
    """

    hduls_list = [hdul for hdul in hduls]
    outputs = []

    if reference is None:
        reference = snr.snr(hduls_list, name)[name]
    # click.echo(reference.header["ORIGNAME"])
    # FIXME log ref name
    np_ref = reference
    try:
        np_ref = np_ref.data
    except AttributeError:
        pass

    for i,hdul in enumerate(hduls_list): ## iterating through list of hdul's. Need to find SCI HDU in stack. 
        np_src = hdul[name] # np_src = hdul["SCI"]
         
        # possibly unneccessary but unsure about scoping
        output = np.array([])
        
        try:
            output = astroalign.register(np_src, np_ref)[0]
        except:
            np_src = hdul[name].data.byteswap().newbyteorder()
            output = astroalign.register(np_src, np_ref)[0]
            pass

        if hasattr(hdul[name], "data"): ##not sure what happens if this if statement fails. previously it would just write empty data.
            idx = hdul.index_of(name)
            hdul[idx].data = output
            hdul[idx].header['EXTNAME'] = ("ALGN    ")

    return (hdul for hdul in hduls_list)

@cli.cli.command("align")
@click.option("-n", "--name", default="SCI", help="The HDU to be aligned.")
@cli.operator
#TODO use CAT sources if they exist

## align function wrapper
def align_cmd(hduls, name="SCI", reference=None):
    """
    Aligns the source astronomical image(s) to the reference astronomical image
    \b
    :param hduls: list of fitsfiles
    :return: list of fistfiles with <name> HDU aligned
    """
    return align([hduls for hduls in hduls], name, reference)
