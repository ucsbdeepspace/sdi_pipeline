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
from photutils.aperture import EllipticalAperture
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
    ds = []
    for c in coord:
        d = (c.ra - cone_center.ra) ** 2 + (c.dec - cone_center.dec) ** 2
        # The 0.0001 so we don't get edge effects
        ds.append(d)
    return [d < (cone_radius ** 2)for d in ds]

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
def norm(ref_table):
    gaia_mag_err = 0.16497
    try:
        #gaia reference data
        ref_mag = np.copy(ref_table['phot_g_mean_mag'])
        g_median = np.max(ref_mag[ref_mag >= 0])
        ref_mag[ref_mag == 0] = g_median
        
        g_rp_g = np.copy(ref_table["phot_rp_mean_mag"])
        rp_median = np.max(g_rp_g[g_rp_g >= 0])
        g_rp_g[g_rp_g == 0] = rp_median
        
        g_bp_g = np.copy(ref_table["phot_bp_mean_mag"])
        bp_median = np.max(g_bp_g[g_bp_g >= 0])
        g_bp_g[g_bp_g == 0] = bp_median
        
        #Making sure that all three values are present for calculating transformed mag
        gaia_mag_transformed = []
        for i,j,k in zip(range(0,len(ref_mag)),range(0,len(g_rp_g)),range(0,len(g_bp_g))):
            if 25>ref_mag[i]>0.5  and 25>g_rp_g[j]>0.5  and 25>g_bp_g[k]>0.5:    
                #transforming reference magnitudes to SDSS12 g
                gaia_mag_transformed.append(ref_mag[i]-0.13518+0.4625*(g_bp_g[k]-g_rp_g[j])+0.25171*(g_bp_g[k]-g_rp_g[j])**2-0.021349*(g_bp_g[k]-g_rp_g[j])**3)
            else:
                gaia_mag_transformed.append(np.nan)
        
        gaia_mag_transformed = np.array(gaia_mag_transformed)
    except ValueError:
        gaia_mag_transformed = []
        for table in ref_table:
            #gaia reference data
            ref_mag = table['phot_g_mean_mag']
            g_median = np.max(ref_mag[ref_mag >= 0]) #If zero, it will just fill it in with a zero
            ref_mag[ref_mag == 0] = g_median
            
            g_rp_g = table["phot_rp_mean_mag"]
            rp_median = np.max(g_rp_g[g_rp_g >= 0])
            g_rp_g[g_rp_g == 0] = rp_median
            
            g_bp_g = table["phot_bp_mean_mag"]
            bp_median = np.max(g_bp_g[g_bp_g >= 0])
            g_bp_g[g_bp_g == 0] = bp_median
            
            one_star = []
            #Making sure that all three values are present for calculating transformed mag
            for i,j,k in zip(range(0,len(ref_mag)),range(0,len(g_rp_g)),range(0,len(g_bp_g))):
                if 25>ref_mag[i]>0.5  and 25>g_rp_g[j]>0.5  and 25>g_bp_g[k]>0.5:
                    #transforming reference magnitudes to SDSS12 g
                    one_star.append(ref_mag[i]-0.13518+0.4625*(g_bp_g[k]-g_rp_g[j])+0.25171*(g_bp_g[k]-g_rp_g[j])**2-0.021349*(g_bp_g[k]-g_rp_g[j])**3)
            gaia_mag_transformed.append(one_star)
    
    return gaia_mag_transformed

def photometry(x,y,a,b,sci_img):
    
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
    source_ap = EllipticalAperture(source_pos, a, b)
    source_flux = aperture_photometry(data, source_ap, error=errmap)
    mag = -2.5 * np.log10(source_flux['aperture_sum'] * gain / exp_time)
    magerr = np.sqrt(((-5 / ((2 * np.log(10)) * (source_flux['aperture_sum'] * gain))) * (source_flux['aperture_sum_err'] * gain)) ** 2)
    
    return mag, magerr

#%% 3
#Start by retrieveing radec for all reference stars
files = glob('/home/pkotta/sdi_output_test/*.fits',recursive=True)
ims = [fits.open(f) for f in files]
del files

#Science images:
sci_f = glob('/home/pkotta/GTAnd_SCI/*.fits',recursive=True)
sci_ims = [fits.open(f) for f in sci_f]
del sci_f
#%% 4
try:
    ts = Time([sci_ims[i][0].header['DATE-OBS'] for i in range(len(ims))]).mjd
except IndexError:
    ts = np.array([i['SCI'].header['OBSMJD'] for i in sci_ims])
#%% 5
#Find the comp_stars in the ref hdu
ref_source = hdultocluster(ims, name="REF", tablename="ROBJ")
cat_source = hdultocluster(ims, name = 'CAT', tablename= 'OBJ')
var_source = hdultocluster(ims, name = 'XRT', tablename= 'XOBJ')
#%% 6
ref_ra = ims[0]['REF'].data['ra']
ref_dec = ims[0]['REF'].data['dec']
refcoord = SkyCoord(ref_ra,ref_dec,frame = 'icrs',unit='degree')
del ref_ra, ref_dec
#%% 7
#Then use _in_cone to see which stars are close to the target star
#target_coord =[SkyCoord(11.291,41.508, frame = 'icrs', unit = 'degree'),SkyCoord(11.291,41.508, frame = 'icrs', unit = 'degree'),SkyCoord(11.291,41.508, frame = 'icrs', unit = 'degree'),SkyCoord(11.291,41.508, frame = 'icrs', unit = 'degree'),SkyCoord(11.291,41.508, frame = 'icrs', unit = 'degree'),SkyCoord(11.291,41.508, frame = 'icrs', unit = 'degree')]
cat_ra = ims[0]['CAT'].data['ra'][0:6]
cat_dec = ims[0]['CAT'].data['dec'][0:6]
#print(cat_ra, cat_dec, "viola")
target_coord = SkyCoord(cat_ra,cat_dec, frame = 'icrs', unit = 'degree')
del cat_ra, cat_dec, ims
#%% 8
#Remove target stars from reference coord list
comp_coords = []
for coord in refcoord:
    if np.array(_in_cone(target_coord,coord, u.Quantity(0.002, u.deg))).any() ==  True:
        pass
    else:
        comp_coords.append(coord)
#%% 9
#Find the comp stars in every image

#These are the comp stars from the ref hdus

#Check that the comp_coords are not present in XRT (non-variable sources)
var_sources = [cluster_search_radec(var_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords]
var_coords = [SkyCoord(v['ra'],v['dec'], frame = 'icrs', unit = 'degree') for v in var_sources][0]
del var_sources, var_source
#%% 10
nonvar_idx = []
for v in var_coords:
    for c in range(0,len(comp_coords)):
        if v!=comp_coords[c]:
            nonvar_idx.append(c)
        else:
            pass
comp_coords_new = np.array(comp_coords)[list(set(nonvar_idx))]
del var_coords, nonvar_idx, refcoord, comp_coords, c
#%% 11
ref_in_cat = [cluster_search_radec(cat_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords_new]
ref_in_ref = [cluster_search_radec(ref_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords_new]
targets = [cluster_search_radec(cat_source, t.ra.deg, t.dec.deg) for t in target_coord]
del comp_coords_new
#%% 13
#Now check for zero x, y, and a in the comp stars. If there are zero values, remove the comp source entirely
c_idx = []
for i in range(0,len(ref_in_cat)):
    if ref_in_cat[i]['x'].all()!=0 and ref_in_cat[i]['y'].all()!=0 and ref_in_cat[i]['a'].all()!=0 and ref_in_cat[i]['b'].all()!=0:
        c_idx.append(i)

comp_sources = np.array(ref_in_cat)[c_idx]
comp_stars = np.array(ref_in_ref)[c_idx]
del c_idx, cat_source, ref_source, ref_in_cat, ref_in_ref
#%%
ref_mag = [] #magnitudes for the reference stars in first image
bright_idx = []
mag_temp = norm(comp_stars)
for mag in mag_temp:
    if np.mean(mag) < 19:
        ref_mag.append(mag)
        bright_idx.append(mag_temp.index(mag))
    else:
        pass
ref_mag = [r[0] for r in ref_mag]
comp_sources_new = np.array(comp_sources)[bright_idx]
del bright_idx, mag_temp, comp_sources,  comp_stars, mag
"""For one target star
#Do photometry on all of the ref sources
import timeit
start_time = timeit.default_timer()
target_mag = []
target_magerr = []
for im in sci_ims:
    idx_im = sci_ims.index(im)
    instrumental_mag = []
    ref_in_magerr = []
    for ref in comp_sources_new:
        i_mag = photometry(ref['x'][idx_im],ref['y'][idx_im],ref['a'][idx_im],ref['b'][idx_im],im)
        instrumental_mag.append(i_mag[0][0])
        ref_in_magerr.append(i_mag[1][0])
        del i_mag
    del ref_in_magerr
    idx = np.where([np.isnan(r)==False for r in ref_mag])[0]
    x = [instrumental_mag[i] for i in idx]
    y = [ref_mag[i] for i in idx]
    fit, sum_sq_resid, rank, singular_values, rcond = np.polyfit(x[1:], y[1:], 1, full=True)
    fit_fn = np.poly1d(fit)
    residuals = fit_fn(x)-y
    del rank, singular_values, rcond, sum_sq_resid,x ,y,idx,instrumental_mag
    #fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    #plt.plot(x,y, 'yo', x, fit_fn(x), '--k')
    #plt.show()        
    # fit_fn is now a function which takes in x and returns an estimate for y, 
    #Use the fit from above to calculate the target magnitude
    target_inst_mag = photometry(target['x'][idx_im],target['y'][idx_im], target['a'][idx_im],target['b'][idx_im], sci_ims[idx_im])[0][0]
    target_mag.append(fit_fn(target_inst_mag))
    target_magerr.append(np.std(residuals))
    del im
elapsed = timeit.default_timer() - start_time
"""
#%%
""" For multiple Target stars """
#Do photometry on all of the ref sources
import timeit
start_time = timeit.default_timer()
fits = []
for im in sci_ims:
    idx_im = sci_ims.index(im)
    instrumental_mag = []
    ref_in_magerr = []
    for ref in comp_sources_new:
        i_mag = photometry(ref['x'][idx_im],ref['y'][idx_im],ref['a'][idx_im],ref['b'][idx_im],im)
        instrumental_mag.append(i_mag[0][0])
        ref_in_magerr.append(i_mag[1][0])
        del i_mag
    #commenting out for error del ref_in_magerr
    idx = np.where([np.isnan(r)==False for r in ref_mag])[0]
    x = [instrumental_mag[i] for i in idx]
    y = [ref_mag[i] for i in idx]
    fit, sum_sq_resid, rank, singular_values, rcond = np.polyfit(x[1:], y[1:], 1, full=True)
    fit_fn = np.poly1d(fit)

    #defining the parameter error for the linear fit
    coeff, cov = np.polyfit(x[1:], y[1:],1,cov = 'true')
    parameter_err = np.array(np.sqrt(np.diag(cov)))
    fit_err =np.array(parameter_err[0])
    fit_slope =np.abs(coeff[0])
    fit_intercept =np.abs(coeff[1])
    int_err = np.array(parameter_err[1])

    residuals = fit_fn(x)-y
    del rank, singular_values, rcond, sum_sq_resid,x ,y,idx,instrumental_mag, im
    fits.append(fit)
    #fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    #plt.plot(x,y, 'yo', x, fit_fn(x), '--k')
    #plt.show()
    #fit_fn is now a function which takes in x and returns an estimate for y, 
    #Use the fit from above to calculate the target magnitude
elapsed = timeit.default_timer() - start_time
print(elapsed)
#%%
target_magerr = []
#Applying the fits to each target star
import csv
for idx, target in enumerate(targets):
    t_idx = np.where(target['x']!=0)
    targets_new = np.array(target[t_idx])
    sci_ims_new = [sci_ims[t] for t in t_idx[0]]
    fits_new = [fits[t] for t in t_idx[0]]
    ts_new = ts[t_idx]
    print(targets_new['ra'][0],targets_new['dec'][0])
    target_mag = []
    for im in sci_ims_new:
        idx_im = sci_ims_new.index(im)
        target_inst_mag = photometry(targets_new['x'][idx_im],targets_new['y'][idx_im], targets_new['a'][idx_im],targets_new['b'][idx_im], sci_ims_new[idx_im])[0][0]
        fit = fits_new[idx_im]
        fit_fn = np.poly1d(fit)
        target_mag.append(fit_fn(target_inst_mag))
        #del target_inst_mag
    #error propagation
    #target_inst_mag = photometry(targets_new['x'][idx_im],targets_new['y'][idx_im], targets_new['a'][idx_im],targets_new['b'][idx_im], sci_ims_new[idx_im])[0][0]
        def relative_u(sigma,x):
            return sigma/(np.abs(x))
        slope_u = relative_u(fit_err, fit_slope)
        mag_u = np.mean(relative_u(ref_in_magerr, target_inst_mag))
    #for i in range(0,len(mag_u)):
     #   a = np.sqrt(mag_u[i]**2+slope_u**2)
        def err(mag, slope):
            return np.sqrt(mag**2+slope**2)
        mx_err = err(mag_u, slope_u)
        def corr(err, slope, mag):
            return err*slope*mag
        part_err = corr(mx_err, fit_slope, target_inst_mag)
        def add(err):
            b = int_err
            return np.sqrt(err**2+b**2)/2
        err = add(part_err)
        target_magerr.append(err)
        del target_inst_mag
    #target_magerr.append(err)

    rows = zip(target_mag,target_magerr, ts_new)
    #rows = zip(target_mag, ts_new)
    with open(str(idx)+'.csv', "w") as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)
    print('done!')
