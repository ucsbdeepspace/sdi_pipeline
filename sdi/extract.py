"""
Extracts sources from a list of hduls.
History:
    Created by 2019-09-05
        Varun Iyer <varun_iyer@ucsb.edu>
    Refactored for 2.0.0 on 2021-04-23
        Varun Iyer <varun_iyer@ucsb.edu>
"""
import click
import sep
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from . import _cli as cli
import time 


def extract(hduls, stddev_thresh=3.0, read_ext="SUB", write_ext="XRT"):
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
    print(" ")
    print("Beginning Extraction")
    start = time.perf_counter()
    for hdul in hduls:
        data = hdul[read_ext].data
        bkg = None
        try:
            bkg = sep.Background(data)
        except ValueError:
            bkg = sep.Background(data.byteswap().newbyteorder())
        sources = sep.extract(data - bkg.back(), bkg.globalrms * stddev_thresh,
                              segmentation_map=False)
        # Convert xy to RA/DEC in ICRS
        coords = wcs.utils.pixel_to_skycoord(sources['x'], sources['y'], wcs.WCS(hdul[read_ext].header), origin=0)
        coords_icrs = SkyCoord(coords, frame='icrs')
        ra = coords_icrs.ra.degree
        dec = coords_icrs.dec.degree
        extname = write_ext
        extver = None
        header = fits.Header([fits.Card("HISTORY", "Extracted by sep.")])
        try:
            extver = int(write_ext[1])
            extname = write_ext[0]
        except (ValueError, TypeError):
            pass
        cat = fits.BinTableHDU(data=sources, header=header,
                               name=extname, ver=extver)
        print('name=', extname)
        out_cols = cat.data.columns
        new_cols = fits.ColDefs([fits.Column(name='ra', format='D', array=ra),
                                 fits.Column(name='dec', format='D', array=dec)])
        new_table = fits.BinTableHDU(header=header, name=extname, ver=extver).from_columns(out_cols+new_cols)
        new_hdu = fits.BinTableHDU(new_table.data, header=header, name=extname, ver=extver)
        hdul.append(new_hdu)
        stop = time.perf_counter()
        print("Extraction Complete")
        print("Time Elapsed: {} sec".format(round(stop-start, 4)))
        yield hdul


@cli.cli.command("extract")
@click.option("-t", "--threshold", default=3.0,
              help="A threshold value to use for source extraction in terms of"
              "the number of stddevs above the background noise.", type=float)
@click.option("-r", "--read_ext", default="SUB",
              help="An index number or ext name that identifies the data in"
              "input hduls that you want source extraction for. For LCO, this "
              "is 1 or SUB.")  # We might want to re-word this now.
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
