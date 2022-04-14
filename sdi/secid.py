"""
Extracts sources from a list of hduls.
History:
    Created by Alex Thomas <agthomas@ucsb.edu>
        2021-09-13
"""
import click
from astropy.coordinates import SkyCoord
from astropy import units as u
from . import _cli as cli

RA_BOUNDS = [12.44, 11.97, 11.48, 11, 10.52, 10.04, 9.56]   #Bounds used to define sections
_RA_SCALE = [-1, 0, 9, 18, 27, 36, 45]  #Values to sum to find section ID according to a grid
DEC_BOUNDS = [42.77, 42.45, 42.13, 41.81, 41.49, 41.17, 40.85, 40.53, 40.21, 39.89]
_DEC_SCALE = [-1, 9, 8, 7, 6, 5, 4, 3, 2, 1]
def secid(hduls, read_ext='SCI', write_ext=-1):
    """
    secid adds M31 section ID information to the header of the HDU specified in write_ext.
    If an HDUL is not in M31, the SECID is set to -1.
    :param hduls: a collection or generator of HDUL to process.
    :param read_ext: an extname of an HDU to read ra and dec information from.
    :param write_ext: an extname or index number of an HDU to write secid information into.
    """
    if write_ext == -1:     #Defaults the write_ext to the read_ext
        write_ext = read_ext
    for hdul in hduls:
        RA = hdul[read_ext].header['RA']    #Reads in RA and dec
        DEC = hdul[read_ext].header['DEC']
        sky_coord = SkyCoord(RA+' '+DEC, uni=(u.hourangle, u.deg))
        sc_ra = sky_coord.ra*u.deg  #Combines RA and Dec into a skycoord with correct units
        sc_dec = sky_coord.dec*u.deg
        ra_id = -1   #starts the values at -1 in case they are too small
        dec_id = -1
        #Constraints the SECID to one 'row' based on the ra
        for i in enumerate(RA_BOUNDS):
            if sc_ra.value > RA_BOUNDS[i]: #Compares value to bounds to find section.
                ra_id = _RA_SCALE[i]
                break

        #Figures out what SECID you are precisely in using dec
        for j in enumerate(DEC_BOUNDS):
            if sc_dec.value > DEC_BOUNDS[j]:
                dec_id = _DEC_SCALE[j]
                break

        if dec_id != -1 and ra_id != -1: #Creates the sec_id
            sec_id = dec_id + ra_id
        else:
            sec_id = -1

        #Writes the secid into the hdul
        hdul[write_ext].header.set('SECID', sec_id, comment='Setion ID in M31 as according to Lubin mosaic')
        yield hdul  #Passes the hdul down the pipeline

@cli.cli.command("secid")
@click.option("-r", "--read_ext", default='SCI',
              help="An index number or ext name that identifies the data in"
              "input hduls that you want source extraction for. For LCO, this "
              "is 0 or SCI.")
@click.option("-w", "--write_ext", default=-1,
              help="An extension name and extension version that will identify"
              "the HDUL that the resulting secid gets written to. Default"
              "is -1")
@cli.operator

def secid_cmd(hduls, read_ext, write_ext):
    """
    Use secid to find the section of an image
    """
    return secid(hduls, read_ext, write_ext)
