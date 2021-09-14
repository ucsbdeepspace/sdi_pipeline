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

RABOUND = [12.44,11.97,11.48,11,10.52,10.04,9.56]    
_RASCALE = [-1,0,9,18,27,36,45]
DECBOUND = [42.77,42.45,42.13,41.81,41.49,41.17,40.85,40.53,40.21,39.89]
_DECSCALE = [-1,9,8,7,6,5,4,3,2,1]
def secid(hduls, read_ext='SCI', write_ext=-1):
    """
    secid adds M31 section ID information to the header of the HDU specified in write_ext.
    If an HDUL is not in M31, the SECID is set to -1.
    :param hduls: a collection or generator of HDUL to process.
    :param read_ext: an extname of an HDU to read RA and DEC information from.
    :param write_ext: an extname or index number of an HDU to write secid information into.
    """
    if write_ext == -1:
        write_ext = read_ext
    for hdul in hduls:
        #print(hduls)
        #print(hdu)
        RA = hdul[read_ext].header['RA']
        DEC = hdul[read_ext].header['DEC'] 
        sc = SkyCoord(RA + ' ' + DEC, unit = (u.hourangle,u.deg))
        SCRA = sc.ra*u.deg
        SCDEC = sc.dec*u.deg
        RAid = -1
        DECid = -1
        #Constraints the SECID to one 'row' based on the RA
        for i,j in enumerate(RABOUND):
            if SCRA.value > RABOUND[i]:
                RAid = _RASCALE[i]
                break
        
        #Figures out what SECID you are precisely in using DEC
        for i,j in enumerate(DECBOUND):
            if SCDEC.value > DECBOUND[i]:
                DECid = _DECSCALE[i]
                break
        #Creates the SECid if both the RA and DEC were found within bounds.
        if DECid != -1 and RAid != -1: 
            SECid = DECid + RAid
        else:
            SECid = -1
        #print(SECid)
        hdul[write_ext].header.set('SECID', SECid)
        #print(hduls[write_ext].header['SECID'])
        yield hdul

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
