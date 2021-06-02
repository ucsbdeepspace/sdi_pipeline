from astropy.io import fits
import numpy as np
import scipy.optimize as scp

#Function for calculating apparent magnitudes
def normalize(flux,a,b):
    return a*np.log10(flux)+b

def norm_sources(hduls, source_ext='CAT', ref_ext='REF', sel_method: ("all", "nrandom", "interactive","test","gaia") = "all", n=None):
    """
    This function appends an 'refmag' and 'refmagerr' field to a catalog.
    These fields contain a magnitude derived by comparing sources with reference stars 
    from remote databases and calculating a relationship between fluxes in each image and 
    magnitudes from reference stars,
    :param hduls: a collection of HDUL to normalize
    :param source_ext: an extname or index of the HDU in HDUL to use as raw source data
    :param ref_ext: an extname or index which specifies an HDU of remote reference star data
    :param sel_method: how to choose the reference stars to use. options are all, nrandom, or
        interactive (which lands you in an interpreter where you can manually select stars)
    :param n: number of refstars to use (only used with "nrandom" option
    :returns: HDULs with 'refmag' and 'refmagerr' fields in the source_ext HDU
    """
    #### NEED TO PUT THE WHOLE THING IN A FOR LOOP TO DEAL WITH EACH FILE
    #Input hduls are the numbered hduls after running ref on the images.
    cat_table = fits.open(hduls)[source_ext]
    ref_table = fits.open(hduls)[ref_ext]
    
    #Extracting flux and mag data from the reference and cat tables within the HDUL
    cat_flux = cat_table.data['flux']
    try:
        ref_flux = ref_table.data['phot_g_mean_flux']
    except:
        raise ValueError('Error fetching flux from REF file: No phot_g_mean_flux data in the referenec star information')
    try:
        ref_mag = ref_table.data['phot_g_mean_mag']
    except:
        raise ValueError('Error fetching mag from REF file: No phot_g_mean_mag data in the reference star information')
    #Transforming ref_mag and ref_flux to work with PTF data
    """ PTF used SDSS g filter, and g-r filter color """
    """ERROR ANALYSIS NOT YET IMPLEMENTED"""
    #If filter data is available in images
    if sel_method=='gaia':
        g_rp_f = ref_table.data["phot_rp_mean_flux"]
        g_bp_f = ref_table.data["phot_bp_mean_flux"]
        gaia_flux_transformed = ref_flux-0.13518+0.4625*(g_bp_f-g_rp_f)+0.25171*(g_bp_f-g_rp_f)**2-0.021349*(g_bp_f-g_rp_f)**3
        g_rp_g = ref_table.data["phot_rp_mean_mag"]
        g_bp_g = ref_table.data["phot_bp_mean_mag"]
        gaia_mag_transformed = ref_mag-0.13518+0.4625*(g_bp_g-g_rp_g)+0.25171*(g_bp_g-g_rp_g)**2-0.021349*(g_bp_g-g_rp_g)**3
        coeff,pcov =scp.curve_fit(normalize,gaia_flux_transformed,gaia_mag_transformed)
        a,b = coeff
        mag = normalize(cat_flux,a,b)
    #Method for a random number of reference stars
    if sel_method == "nrandom":
        index = np.random.randint(0,len(ref_flux)-1)
        coeff,pcov = scp.curve_fit(normalize,ref_flux[:index],ref_mag[:index])
        a,b = coeff
        mag = normalize(cat_flux,a,b)
    #Method for using all the reference stars
    if sel_method == "all":
        coeff,pcov = scp.curve_fit(normalize,ref_flux,ref_mag)
        a,b = coeff
        mag = normalize(cat_flux,a,b)
    #Method for the interactive
    #if sel_method == "interactive":
    #not implemented
    if sel_method == "test":
        mag = [ref_mag[0]-2.512*(f/ref_flux[0]) for f in cat_flux]
    #Raise value error if none of the above (use elif)
    else:
        raise ValueError('No valid selection method was selected')
    
    # Writing the magnitudes and magerrors to the source_ext file in the HDUL
    
    #First build the output table
    out_cols = cat_table.columns
    new_cols = fits.ColDefs([fits.Column(name='mag', format='D',array=mag),fits.Column(name='magerr', format='D',array=np.zeros(len(cat_table)))])
    hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    hduls.append(hdu)
    return hduls
