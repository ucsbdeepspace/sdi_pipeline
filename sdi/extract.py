"""
Extracts sources from a list of hduls.
History:
    Created by 2019-09-05
        Varun Iyer <varun_iyer@ucsb.edu>
    Refactored for 2.0.0 on 2021-04-23
        Varun Iyer <varun_iyer@ucsb.edu>
"""
import click
from astropy.io import fits
import sep
from . import _cli as cli


def extract(hduls, stddev_thresh=3.0, read_ext=0, write_ext="XRT"):
    """
    Uses sep to find sources on a residual image(s)
    :param hduls: a list of HDUL to use as science data
    :param thresh: a threshold value for source extraction in terms of
        # of stddevs above background noise
    :param read_ext: the HDUL index to use as the base image. Default is 0. For
        LCO, 'SCI' would do the trick; for others 'PRIMARY'.
    :param write_ext: the HDUL index to write the catalog to. Default is 'XRT'
        for eXtRacT (close enough). Can also be a tuple (extname, extver)
    :return: a list of HDUL, each with a new `write_ext` HDU appended that is a
        record of extracted sources
    """
    for hdul in hduls:
        data = hdul[read_ext].data
        bkg = None
        try: 
            bkg = sep.Background(data)
        except:
            data = data.byteswap().newbyteorder()
            bkg = sep.Background(data)
        sources = sep.extract(data - bkg.back(), bkg.globalrms * stddev_thresh,
                       segmentation_map=False)
        extname = write_ext
        extver = None
        header = fits.Header([fits.Card("HISTORY", "Extracted by sep.")])
        try:
            extver = int(write_ext[1])
            extname = write_ext[0]
        except (ValueError, TypeError):
            pass
        hdul.append(fits.BinTableHDU(data=sources, header=header,
                                     name=extname, ver=extver))
        yield hdul

@cli.cli.command("extract")
@click.option("-t", "--threshold", default=3.0,
              help="A threshold value to use for source extraction in terms of"
              "the number of stddevs above the background noise.", type=float)
@click.option("-r", "--read_ext", default=0,
              help="An index number or ext name that identifies the data in"
              "input hduls that you want source extraction for. For LCO, this "
              "is 0 or SCI.")
@click.option("-w", "--write_ext", default=("XRT", 1),
              help="An extension name and extension version that will identify"
              "the HDUL that the resulting BinTable gets written to. Default"
              "is `XRT 1`", type=(str, int))

@cli.operator
def extract_cmd(hduls, threshold, read_ext, write_ext):
    """
    Uses sep to find sources in ImageHDU data.
    """
    try:
        read_ext = int(read_ext)
    except ValueError:
        pass
    return extract(hduls, threshold, read_ext, write_ext)
