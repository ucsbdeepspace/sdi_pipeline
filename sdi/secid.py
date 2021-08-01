from astropy.coordinates import SkyCoord
from astropy import units as u

RABOUND = [12.44,11.97,11.48,11,10.52,10.04,9.56]    
_RASCALE = [-1,0,9,18,27,36,45]
DECBOUND = [42.77,42.45,42.13,41.81,41.49,41.17,40.85,40.53,40.21,39.89]
_DECSCALE = [-1,9,8,7,6,5,4,3,2,1]
def secid(hduls, read_ext='SCI', write_ext=-1):
    """
    secid adds M31 section ID information to the header of the HDU specified in write_ext.
    If an HDUL does not have WCS data or is not in M31, a ValueError will be raised.
    :param hduls: a collection or generator of HDUL to process.
    :param read_ext: an extname or index number of an HDU to pull RA and DEC information from.
    :param write_ext: an extname or index number of an HDU to write secid information into.
    """
    if write_ext == -1:
        write_ext = read_ext
    #print(hduls)
    primary = hduls[read_ext]
    #print(primary)
    RA = primary.header['RA']                                                           
    DEC = primary.header['DEC']                                                         
    sc = SkyCoord(RA + ' ' + DEC, unit = (u.hourangle,u.deg))                       
    SCRA = sc.ra*u.deg                                                             
    SCDEC = sc.dec*u.deg   
    RAid = -1
    DECid = -1
    #Constraints the SECID to one 'row' based on the RA
    for i in range(len(RABOUND)):
        if SCRA.value > RABOUND[i]:
            RAid = _RASCALE[i]
            break
    
    #Figures out what SECID you are precisely in using DEC
    for i in range(len(DECBOUND)):
        if SCDEC.value > DECBOUND[i]:
            DECid = _DECSCALE[i]
            break
    #Creates the SECid if both the RA and DEC were found within bounds.
    if DECid != -1 and RAid != -1: 
        SECid = DECid + RAid
    else:
        SECid = -1
    #print(SECid)
# if HDUL does not have WCS data OR SECid == -1:
#    raise ValueError("bad SECid")
## should this be one error or one error in each case?
##What is WCS data?
    hduls[write_ext].header.append(('SECID', SECid))
    #print(hduls[write_ext].header['SECID'])

    yield hduls
