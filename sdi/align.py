"""
align -- this module aligns sets of astronomical data
HISTORY
    Created/Substantially refactored on 2019-09-01
        Varun Iyer <varun_iyer@ucsb.edu>
"""
# general imports

import click
from . import _cli as cli
import numpy as np
from astropy.io.fits import PrimaryHDU, HDUList
#from sources import Source
import astroalign
from _scripts import snr

def align(hduls, name="SCI", reference=None):
    """ 
    Aligns the source astronomical image(s) to the reference astronomical image 
    \b 
    :param hduls: list of fitsfiles 
    :return: list of fistfiles with <name> HDU aligned 
    """ 

    hduls_list = [hdul for hdul in hduls]
    sources = [hdul[name] for hdul in hduls_list]
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

    for source in sources:
        np_src = source
        try:
            np_src = source.data
        except AttributeError:
            pass
        # possibly unneccessary but unsure about scoping
        output = np.array([])

        output = astroalign.register(np_src, np_ref)[0]
        if hasattr(source, "data"):
            output = PrimaryHDU(output, source.header)
        outputs.append(HDUList([output]))

    return (hdul for hdul in outputs)

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
    return align([hduls for hduls in hduls],name,reference)
