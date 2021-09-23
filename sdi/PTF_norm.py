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
from scipy.optimize import curve_fit

def normalize(flux,ref_flux,ref_mag):
    #If ref flux is longer
    try:
        return ref_mag[0:len(flux)]-2.5*np.log10(flux/ref_flux[0:len(flux)])
    #If cat flux is longer
    except ValueError:
        return ref_mag[0:len(ref_flux)]-2.5*np.log10(flux[0:len(ref_flux)]/ref_flux)

"""
This function takes in one fits file at a time, and normalizes it with respect to the Gaia reference table magnitudes of stars after transforming the Gaia data color information to the SDSS12 g-filter according to the Gaia documentation available at: https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
"""
#TODO: find a way to write information back to the CAT table in the hduls.

def norm(hduls):
    cat_table = hduls['CAT']
    ref_table = hduls['REF']
    pixel_value = hduls['SCI'].data
    gain = hduls['SCI'].header['GAIN']
    try:
        noise = hduls['SCI'].header['FBIAS']
    except:
        noise = hduls['SCI'].header['RDNOISE']
    #Cat data
    cat_flux = cat_table.data['flux']

    #Using the definition from SExtractor documentation
    cat_fluxerr=[]
    xmin = cat_table.data['xmin']
    xmax = cat_table.data['xmax']
    ymin = cat_table.data['ymin']
    ymax = cat_table.data['ymax']
    for a,b,c,d in zip(xmin,xmax,ymin,ymax):
        pixel_object=np.hstack((pixel_value[c][a:b],pixel_value[d][a:b]))
        cat_fluxerr.append(np.sqrt(np.sum(noise**2+pixel_object/gain)))
    cat_fluxerr = np.array(cat_fluxerr)
    #gaia reference data
    ref_mag = ref_table.data['phot_g_mean_mag']
    g_rp_g = ref_table.data["phot_rp_mean_mag"]
    g_bp_g = ref_table.data["phot_bp_mean_mag"]
    
    #Making sure that all three values are present for calculating transformed mag
    gaia_mag_transformed = []
    flux_temp = [] #is either equal or smaller than gaia_mag_transformed
    indices = []
    for i,j,k in zip(range(0,len(ref_mag)-1),range(0,len(g_rp_g)-1),range(0,len(g_bp_g)-1)):
        if 25>ref_mag[i]>0.5  and 25>g_rp_g[j]>0.5  and 25>g_bp_g[k]>0.5:    
            #transforming reference magnitudes to SDSS12 g
            gaia_mag_transformed.append(ref_mag[i]-0.13518+0.4625*(g_bp_g[k]-g_rp_g[j])+0.25171*(g_bp_g[k]-g_rp_g[j])**2-0.021349*(g_bp_g[k]-g_rp_g[j])**3)
            try:
                flux_temp.append(cat_flux[i])
                indices.append(i)
            except IndexError:
                pass
    gaia_mag_err = 0.16497
    
    gaia_mag_transformed = np.array(gaia_mag_transformed)
    flux_temp = np.array(flux_temp)
    #transforming reference flux to SDSS12 g    
    #Find c by fitting to our data. It is linear in c and not a very complicated model so we can use np.linalg.lstsq
    #If gaia is shorter
    try:
        A = np.ones((len(gaia_mag_transformed),1))
        A[:, 0] = 10**(gaia_mag_transformed/(-2.5))
        B = []
        for i in A:
            if i != np.inf:
              B.append(i)
        B = np.array(B)
        c_best = np.linalg.lstsq(B, flux_temp[0:len(B)], rcond=None)[0]
    #If cat_flux is shorter
    except:
        A = np.ones((len(flux_temp),1))
        A[:, 0] = 10**(gaia_mag_transformed[0:len(flux_temp)]/(-2.5))
        B = []
        for i in A:
            if i != np.inf:
                B.append(i)
        B = np.array(B)
        c_best = np.linalg.lstsq(B, flux_temp[0:len(B)], rcond=None)[0]
        

    gaia_flux_transformed = np.dot(B, c_best)
    c_err = np.sqrt(np.sum(10**(gaia_mag_transformed/(-2.5))*gaia_mag_err))
    #If cat_flux is shorter
    
    cat_mag = np.ones((len(cat_flux),1))
    #Now, we have to find the best fit coefficients for the reference magnitudes
    mag_temp = normalize(flux_temp,gaia_flux_transformed,gaia_mag_transformed)
    cat_mag[np.ix_(indices),0] = mag_temp
    not_indices=(np.where(cat_mag==1)[0])
    #import pdb;pdb.set_trace()
    #Once we get cat_mag, maybe plot the fluxes vs. mags to find the appropriate values of the flux.
    def func(x, b):
        return -2.512*np.log10(x) + b
    coeff, pcov = curve_fit(func, flux_temp, mag_temp)
    cat_mag[np.ix_(not_indices),0] = func(cat_flux[np.ix_(not_indices)],coeff[0])
    perr = np.sqrt(np.diag(pcov))
    #Error analysis: We want to find the error in the magnitude we just found.
    #The flux, gaia-flux and gaia-magnitudes will all contribute to the error in the magnitude.
    
    #First we have to find the ref_flux error
    gaia_flux_err = np.sqrt((10**(gaia_mag_transformed/(-2.5))*c_err)**2+(-(np.log(10)*c_best*gaia_mag_err)/(2.5*10**(gaia_mag_transformed)*(-0.4)))**2)
    
    #Uncertainty in calculated magnitudes:
    #If ref flux is longer
    mag_err = np.ones((len(cat_flux),1))
    import pdb; pdb.set_trace()
    try:
        mag_err[np.ix_(indices),0]=np.sqrt((gaia_mag_err)**2+((-2.5*cat_fluxerr[np.ix_(indices)])/(np.log(10)*flux_temp))**2+((2.5*gaia_flux_err[0:len(flux_temp)])/(np.log(10)*gaia_flux_transformed[0:len(flux_temp)]))**2)
    #If cat flux is longer
    except ValueError:
        mag_err[np.ix_(indices),0]=np.sqrt((gaia_mag_err)**2+((-2.5*cat_fluxerr[np.ix_(indices)])/(np.log(10)*flux_temp[0:len(gaia_flux_transformed)]))**2+((2.5*gaia_flux_err[0:len(gaia_flux_transformed)])/(np.log(10)*gaia_flux_transformed))**2)
    
    mag_err[np.ix_(not_indices),0]=np.sqrt((-2.512/(np.log(10)*cat_flux[np.ix_(not_indices)])*cat_fluxerr[np.ix_(not_indices)])**2+(perr)**2)
    #write cat_mag and mag_err to the hdul
    #Have to redo writing method so it writes to CAT 
    
    #import pdb; pdb.set_trace()
    out_cols = cat_table.columns
    new_cols = fits.ColDefs([fits.Column(name='mag', format='D',array=cat_mag),fits.Column(name='magerr', format='D',array=mag_err),fits.Column(name='fluxerr',format='D',array=cat_fluxerr)])
    hdu = fits.BinTableHDU().from_columns(out_cols + new_cols)
    hduls.append(hdu)
    return hduls
