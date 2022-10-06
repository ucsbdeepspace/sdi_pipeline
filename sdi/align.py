"""
align -- this module aligns sets of astronomical data
HISTORY
    Created/Substantially refactored on 2019-09-01
        Varun Iyer <varun_iyer@ucsb.edu>
    Refactored SNR section on 2022-2-22
        Jasper Webb <jlw@ucsb.edu>
"""
import click
import numpy as np
import astroalign
import time
from _scripts import snr
# from _scripts import snr
# from astropy.io.fits import CompImageHDU
# from sources import Source
from . import _cli as cli

def align(hduls, name="SCI", ref = None):
    """
    Aligns the source astronomical image(s) to the reference astronomical image
    \b
    :param hduls: list of fitsfiles
    :param name: header containing image data.
    :param reference: index of refence image. Defaults to index of image with highest SNR
    :return: list of FITS files with <name> HDU aligned
    """

    print(" ")
    print("Beginning Alignment")
    start = time.perf_counter()
    fits_files = hduls
    hduls_list = [hdul for hdul in hduls]

    # No reference index given. we establish reference based on best signal to noise ratio
    if ref is None:
        print(" ")
        print("Determining Image with greatest SNR")
        #SNR function appends snr to each hdul's "CAT" table and returns list of hduls
        hduls_list = list(snr(hduls_list, name))
        reference = hduls_list[0][name]  # 0th index reference is used by default
        ref_snr = hduls_list[0]["CAT"].header["SNR"]

        for hdul in hduls_list:  # loops though hdul lists finding the hdul with greatest snr
            if hdul["CAT"].header["SNR"] > ref_snr:  # compares SNR value of current hdul to ref
                ref_snr = hdul["CAT"].header["SNR"]
        reference = hdul[name]

    else:  # ref index is provided
        ref = int(ref)
        reference = hduls_list[ref][name]

    try:
        ref_data = reference.data
    except AttributeError:
        print("The reference file have doesn't have Attribute: Data")
    
    # ref_data = hduls_list[0][name].data
    # reference = hduls_list[0][name]
    
    for i, hdul in enumerate(hduls_list):
        np_src = hdul[name].data
        output = np.array([])
        try:
            output = astroalign.register(np_src, ref_data)[0]
        except ValueError:
            np_src = hdul[name].data.byteswap().newbyteorder()
            output = astroalign.register(np_src, ref_data)[0]

        if hasattr(hdul[name], "data"):
            idx = hdul.index_of(name)
            hdul[idx].data = output
            hdul[idx].header['EXTNAME'] = ("ALGN    ")
            hdul[idx].header = reference.header
    stop = time.perf_counter()
    print("Alignment Complete")
    print("Time Elapsed: {} sec".format(round(stop-start, 4)))

    return (hdul for hdul in hduls_list)


@cli.cli.command("align")
@click.option("-n", "--name", default="SCI", help="The HDU to be aligned.")
@click.option("-r", "--ref", default=None, help="Reference Image either the image number or image with highest SNR if None.")
@cli.operator
# TODO: use CAT sources if they exist

# align function wrapper
def align_cmd(hduls, name="SCI", ref=None):
    """
    Aligns the source astronomical image(s) to the reference astronomical image
    \b
    :param hduls: list of fitsfiles
    :return: list of fistfiles with <name> HDU aligned
    """
    return align([hduls for hduls in hduls], name, ref)
