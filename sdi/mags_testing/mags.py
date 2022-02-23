# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 20:52:51 2022

@author: prera
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry
from glob import glob
from photutils.segmentation import make_source_mask
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
#Inter script functions
from collate import hdultocluster, cluster_search_radec

#%%
def _in_cone(coord: SkyCoord, cone_center: SkyCoord, cone_radius: u.degree):
    """
    Checks if SkyCoord coord is in the cone described by conecenter and
    cone_radius
    """
    d = (coord.ra - cone_center.ra) ** 2 + (coord.dec - cone_center.dec) ** 2
    # The 0.0001 so we don't get edge effects
    return d < (cone_radius ** 2)

def find_ref(target: SkyCoord, refcoord: SkyCoord, threshold: u.degree = u.Quantity(0.03, u.deg)):
    ref_comp = []
    for coord in refcoord:
        try:
            for t in target:
                if _in_cone(t, coord, threshold):
                    ref_comp.append(coord)
                else:
                    pass
        except TypeError:
            if _in_cone(target, coord, threshold):
                    ref_comp.append(coord)
            else:
                pass
    
    #This gets us comparison stars, and the star itself
    #Getting rid of the target stars inthe comp list
    thresh = u.Quantity(0.001, u.deg)
    comp = []
    for coord in ref_comp:
        try:
            for t in target:
                if _in_cone(t, coord, thresh):
                    pass
                else:
                    comp.append(coord)
        except TypeError:
            if _in_cone(target, coord, thresh):
                    pass
            else:
                comp.append(coord)
    return comp

#Convert their magnitudes to SDSS g'
def norm(ref_table,cat_table):
    #gaia reference data
    ref_mag = ref_table['phot_g_mean_mag']
    g_rp_g = ref_table["phot_rp_mean_mag"]
    g_bp_g = ref_table["phot_bp_mean_mag"]
    
    gaia_mag_err = 0.16497
    #Making sure that all three values are present for calculating transformed mag
    gaia_mag_transformed = []
    for i,j,k in zip(range(0,len(ref_mag)),range(0,len(g_rp_g)),range(0,len(g_bp_g))):
        if 25>ref_mag[i]>0.5  and 25>g_rp_g[j]>0.5  and 25>g_bp_g[k]>0.5:    
            #transforming reference magnitudes to SDSS12 g
            gaia_mag_transformed.append(ref_mag[i]-0.13518+0.4625*(g_bp_g[k]-g_rp_g[j])+0.25171*(g_bp_g[k]-g_rp_g[j])**2-0.021349*(g_bp_g[k]-g_rp_g[j])**3)
        else:
            gaia_mag_transformed.append(np.nan)
    
    gaia_mag_transformed = np.array(gaia_mag_transformed)
    
    return gaia_mag_transformed

def photometry(x,y,aperture,sci_img):
    
    data = sci_img['SCI'].data
    gain = sci_img['SCI'].header['GAIN']
    exp_time = sci_img['SCI'].header['EXPTIME']

    # Photometry
    # masking
    mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    imgdata_bkgsub = data - median                                       #only use imgdata_bkgsub from here on
    # create an error map (extra steps are taken to handle negative values)
    errmap = np.sqrt(np.sqrt(imgdata_bkgsub**2))/gain
    
    source_pos = np.transpose((x, y))
    source_ap = CircularAperture(source_pos, r=aperture)
    source_flux = aperture_photometry(data, source_ap, error=errmap)
    mag = -2.5 * np.log10(source_flux['aperture_sum'] * gain / exp_time)
    magerr = np.sqrt(((-5 / ((2 * np.log(10)) * (source_flux['aperture_sum'] * gain))) * (source_flux['aperture_sum_err'] * gain)) ** 2)
    
    return mag, magerr

#%%
#Start by retrieveing radec for all reference stars
files = glob('sdi_output_test/*.fits',recursive=True)
ims = [fits.open(f) for f in files]

#Science images:
sci_f = glob('GTAnd_SCI/*.fits',recursive=True)
sci_ims = [fits.open(f) for f in sci_f]
#%%
try:
    t = [sci_ims[i][0].header['DATE-OBS'] for i in range(len(ims))]
    times = Time(t)
    ts = times.mjd
except IndexError:
    t = [i['SCI'].header['OBSMJD'] for i in sci_ims]
    ts = np.array(t)
#%%
ref_ra = ims[0]['REF'].data['ra']
ref_dec = ims[0]['REF'].data['dec']
refcoord = SkyCoord(ref_ra,ref_dec,frame = 'icrs',unit='degree')

#Then use _in_cone to see which stars are close to the target star
target_coord = SkyCoord(11.291,41.508, frame = 'icrs', unit = 'degree')
comp_coords = find_ref(target_coord, refcoord) #in all images
#%%
#Find the comp_stars in the ref hdu
ref_source = hdultocluster(ims, name="REF", tablename="ROBJ")
cat_source = hdultocluster(ims, name = 'CAT', tablename= 'OBJ')
var_source = hdultocluster(ims, name = 'XRT', tablename= 'XOBJ')
#%%
#Find the comp stars in every image

#These are the comp stars from the ref hdus

#Check that the comp_coords are not present in XRT (non-variable sources)
var_sources = [cluster_search_radec(var_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords]
var_coords = [SkyCoord(v['ra'],v['dec'], frame = 'icrs', unit = 'degree') for v in var_sources][0]
#%%
nonvar_idx = []
for v in var_coords:
    for c in range(0,len(comp_coords)):
        if v!=comp_coords[c]:
            nonvar_idx.append(c)
        else:
            pass
comp_coords_new = np.array(comp_coords)[list(set(nonvar_idx))]
#%%
#These are the comp stars in the cat table itself. so these stars are in our images
comp_sources = [cluster_search_radec(cat_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords_new]
comp_stars = [cluster_search_radec(ref_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords_new]
#This is the target star in every image
target = cluster_search_radec(cat_source, target_coord.ra.deg,target_coord.dec.deg)
#%%
#Remove all the images for which the target star could not be found by removing places that have zeroes
t_idx = np.where(target['x']!=0)
target = target[t_idx]
#next remove those indices in all the comp stars and sources
comp_stars = [comp_stars[i][t_idx] for i in range(0,len(comp_stars))]
comp_sources = [comp_sources[i][t_idx] for i in range(0,len(comp_sources))]

#Now check for zero x, y, and a in the comp stars. If there are zero values, remove the comp source entirely
c_idx = []
for i in range(0,len(comp_sources)):
    if comp_sources[i]['x'].all()!=0 and comp_sources[i]['y'].all()!=0 and comp_sources[i]['a'].all()!=0:
        c_idx.append(i)

comp_sources = np.array(comp_sources)[c_idx]
comp_stars = np.array(comp_stars)[c_idx]

#%%
target_mag = []
target_magerr = []
ts = []
#%%
for im in range(0,len(ims)):
    try:
        t = sci_ims[im][0].header['DATE-OBS']
        times = Time(t)
        ts.append(times.mjd)
    except IndexError:
        t = sci_ims[im]['SCI'].header['OBSMJD']
        ts.append(np.array(t))
    #For multiple ref stars
    try:
        #gaia_mag_transformed, gaia_mag_err, gaia_flux_transformed, gaia_flux_err
        ref_mag = [] #magnitudes for the reference stars in first image
        bright_idx = []
        ref_magerr = 0.16497
        for i in range(0,len(comp_stars)):
            mag_temp = norm(comp_stars[i],comp_sources[i])[im] #This can be made better, normalize once outside the for loop and then iterate
            if mag_temp < 18: #This is just for GTAnd. This will change for other stars based on the target star's magnitude
                ref_mag.append(mag_temp)
                bright_idx.append(i)
            else:
                pass
        comp_sources_new = np.array(comp_sources)[bright_idx]
        #Using 'a' as a sloppy alternative to aperture. Maybe look into sep.kron_radius or flux_radius. aperture is a radius.
        #Have to select index [0] for each mag to get just the numerical value without the description of the column object
        
        #For multiple ref stars
        
        #!TODO: Currently, comp_sources still have zeros? np.mean will not work anyway
        instrumental_mag = []
        in_magerr = []
        for i in range(0,len(comp_sources_new)):
            arr = photometry(comp_sources_new[i]['x'][im],comp_sources_new[i]['y'][im], comp_sources_new[i]['a'][im], sci_ims[im])
            instrumental_mag.append(arr[0][0])
            in_magerr.append(arr[1][0])
        
        #Performing the linear fit
        #First find all places where there are non-nan values:
        idx = np.where([np.isnan(r)==False for r in ref_mag])[0]
        x = [instrumental_mag[i] for i in idx]
        y = [ref_mag[i] for i in idx]
        fit, sum_sq_resid, rank, singular_values, rcond = np.polyfit(x[1:], y[1:], 1, full=True)
        fit_fn = np.poly1d(fit)
        residuals = fit_fn(x)-y
        fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
        plt.plot(x,y, 'yo', x, fit_fn(x), '--k')
        plt.show()        
        # fit_fn is now a function which takes in x and returns an estimate for y, 
        #Use the fit from above to calculate the target magnitude
        target_inst_mag = photometry(target['x'][im],target['y'][im], target['a'][im], sci_ims[im])[0][0]
        target_mag.append(fit_fn(target_inst_mag))
        target_magerr.append(np.std(residuals))
    
    except ValueError:
        #!TODO: correct error prop
        #For one ref star, take the star, find the correction, and apply it to the rest of the image
        instrumental_mag = photometry(comp_sources[0]['x'][im],comp_sources[0]['y'][im], comp_sources[0]['a'][im], sci_ims[im])[0][0]
        in_magerr = photometry(comp_sources[0]['x'][im], comp_sources[0]['y'][im], comp_sources[0]['a'][im], sci_ims[im])[1][0]
        corr = np.abs(instrumental_mag)-ref_mag[0]
        target_inst_mag = photometry(target['x'][im],target['y'][im], target['a'][im], sci_ims[im])[0][0]
        target_mag.append(corr + target_inst_mag)
        target_magerr.append(in_magerr)
#%%
rows = zip(target_mag,target_magerr, ts)

import csv

with open('var_bright_test.csv', "w") as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)
