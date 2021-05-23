"""
Convert DEC and RA values to secid for m31
"""
RABOUND = [12.44,11.97,11.48,11,10.52,10.04,9.56]    
_RASCALE = [-1,0,9,18,27,36,45]
DECBOUND = [42.77,42.45,42.13,41.81,41.49,41.17,40.85,40.53,40.21,39.89]
_DECSCALE = [-1,9,8,7,6,5,4,3,2,1]
def secid(hduls, write_ext=0):
    """
    secid adds M31 section ID information to the header of the HDU specified in write_ext.
    If an HDUL does not have WCS data or is not in M31, a ValueError will be raised.
    :param hduls: a collection or generator of HDUL to process.
    :param write_ext: an extname or index number of an HDU to write secid information into.
    """  
      
    RA = sci.header['RA']                                                           
    DEC = sci.header['DEC']                                                         
    sc = SkyCoord(RA + ' ' + DEC, unit = (u.hourangle,u.deg))                       
    SCRA = sc.ra*u.deg                                                              
    SCDEC = sc.dec*u.deg   
    SECid = -1
    RAid = -1
    DECid = -1
    for i in range(len(RABOUND)):
        if RA > RABOUND[i]:
            RAid = _RASCALE[i]
            break
 
    for i in range(len(DECBOUND)):
        if DEC > DECBOUND[i]:
            DECid = _DECSCALE[i]
            break
    if DECid != -1 and RAid != -1: 
        SECid = DECid + RAid

# if HDUL does not have WCS data OR SECid == -1:
#    raise ValueError("bad SECid")
## should this be one error or one error in each case?
##What is WCS data?

    hduls[write_ext].header.set('SECID', SECid)
    print(hduls[write_ext].header['SECID'])
## This doesn't save!

    yield hduls
