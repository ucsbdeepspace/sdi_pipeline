# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 11:52:51 2021

@author: prera

NOTE: Data for this comes from the output fits files of the pipeline, including the REF extension name
"""
from astropy.io import fits
import os
from glob import glob
import numpy as np

def normalize(flux,ref_flux,ref_mag):
    return ref_mag-2.5*np.log10(flux/ref_flux)

"""
This function takes in one fits file at a time, and normalizes it with respect to the Gaia reference table magnitudes of stars after transforming the Gaia data color information to the SDSS12 g-filter according to the Gaia documentation available at: https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
"""
#TODO: find a way to write information back to the XRT table in the hduls. For now it is in CAT

def norm(hduls):
    cat_table = hduls['CAT']
    ref_table = hduls['REF']
    
    #Cat data
    cat_flux = cat_table.data['flux']
    cat_fluxerr = np.std(cat_flux)
    
    #gaia reference data
    ref_mag = ref_table.data['phot_g_mean_mag']
    g_rp_g = ref_table.data["phot_rp_mean_mag"]
    g_bp_g = ref_table.data["phot_bp_mean_mag"]
    
    #transforming reference magnitudes to SDSS12 g
    gaia_mag_transformed = ref_mag-0.13518+0.4625*(g_bp_g-g_rp_g)+0.25171*(g_bp_g-g_rp_g)**2-0.021349*(g_bp_g-g_rp_g)**3
    gaia_mag_err = 0.16497

    #transforming reference flux to SDSS12 g    
    #Find c by fitting to our data. It is linear in c and not a very complicated model so we can use np.linalg.lstsq
    A = np.ones((len(cat_flux),1))
    A[:, 0] = 10**(-2.5*gaia_mag_transformed)
    c_best = np.linalg.lstsq(A, cat_flux, rcond=None)[0]
    gaia_flux_transformed = np.dot(A, c_best)
    
    c_err = np.sqrt(np.sum(10**(2.5*gaia_mag_transformed)*gaia_mag_err))
    
    #Now, we have to find the best fit coefficients for the reference magnitudes
    cat_mag = normalize(cat_flux,gaia_flux_transformed,gaia_mag_transformed)
    
    #Error analysis: We want to find the error in the magnitude we just found.
    #The flux, gaia-flux and gaia-magnitudes will all contribute to the error in the magnitude.
    
    #First we have to find the ref_flux error
    gaia_flux_err = np.sqrt((10**(-2.5*gaia_mag_transformed)*c_err)**2+((np.log(10)*c_best*gaia_mag_err)/(2.5*10**(2.5*gaia_mag_transformed)))**2)
    
    #Uncertainty in calculated magnitudes:
    mag_err = np.sqrt((gaia_mag_err)**2+((-2.5*cat_fluxerr)/(np.log(10)*cat_flux))**2+((2.5*gaia_flux_err)/(np.log(10)*gaia_flux_transformed))**2)

    #write cat_mag and mag_err to the hdul
    out_cols = cat_table.columns
    new_cols = fits.ColDefs([fits.Column(name='mag', format='D',array=cat_mag),fits.Column(name='magerr', format='D',array=mag_err)
    hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    hduls.append(hdu)
    return hduls
