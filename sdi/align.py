"""
align -- this module aligns sets of astronomical data
HISTORY
    Created/Substantially refactored on 2019-09-01
        Varun Iyer <varun_iyer@ucsb.edu>
    Refactored SNR section on 2022-2-22
        Jasper Webb <jlw@ucsb.edu>
"""
# general imports
import click
import numpy as np
from astropy.io.fits import CompImageHDU

# from sources import Source
import astroalign
from .snr_function import snr

# from _scripts import snr
from . import _cli as cli


def align(hduls, name="SCI", ref=None):
    """
    Aligns the source astronomical image(s) to the reference astronomical image
    \b
    :param hduls: list of fitsfiles
    :param name: header containing image data.
    :param reference: index of desired
    :return: list of FITS files with <name> HDU aligned
    """

    hduls_list = [hdul for hdul in hduls]
    outputs = []

    hduls_list = snr(
        hduls_list, name
    )  # SNR function appends snr to each hdul's "CAT" table.
    hduls_list = list(hduls_list)  # convert generator to list.

    if (
        ref is None
    ):  # No reference index given. we establish reference based on best signal to noise ratio

        reference = hduls_list[0][name]  # 0th index reference is used by default
        ref_SNR = hduls_list[0]["CAT"].header["SNR"]

        for (
            hdul
        ) in hduls_list:  # loops though hdul lists finding the hdul with greatest snr
            if (
                hdul["CAT"].header["SNR"] > ref_SNR
            ):  # compares SNR value of current hdul to ref
                ref_SNR = hdul["CAT"].header["SNR"]
        reference = hdul[name]

    else:  # ref index is provided
        reference = hduls_list[ref][name]

    try:
        ref_data = reference.data
    except AttributeError:
        print("The reference file have doesn't have Attribute: Data")
        pass

    for i, hdul in enumerate(hduls_list):
        np_src = hdul[name]

        # possibly unneccessary but unsure about scoping
        output = np.array([])

        try:
            output = astroalign.register(np_src, ref_data)[0]
        except ValueError:
            np_src = hdul[name].data.byteswap().newbyteorder()
            output = astroalign.register(np_src, ref_data)[0]
            pass

        if hasattr(hdul[name], "data"):
            idx = hdul.index_of(name)
            hdul[idx].data = output
            hdul[idx].header["EXTNAME"] = "ALGN    "
            hdul[idx].header = reference.header  # hdul['sci'].header = reference.header

    return (hdul for hdul in hduls_list)


@cli.cli.command("align")
@click.option("-n", "--name", default="SCI", help="The HDU to be aligned.")
@cli.operator
# TODO use CAT sources if they exist

## align function wrapper
def align_cmd(hduls, name="SCI", ref=None):
    """
    Aligns the source astronomical image(s) to the reference astronomical image
    \b
    :param hduls: list of fitsfiles
    :return: list of fistfiles with <name> HDU aligned
    """
    return align([hduls for hduls in hduls], name, ref)
