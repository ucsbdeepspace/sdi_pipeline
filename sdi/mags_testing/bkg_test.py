
"""
Created on Wed Apr 27 2022

@author: Edgar
"""

import numpy as np
import matplotlib.pyplot as plt
import sep
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


mags_arr = []

def _in_cone(coord: SkyCoord, cone_center: SkyCoord, cone_radius: u.degree):
    """
    Checks if SkyCoord coord is in the cone described by conecenter and
    cone_radius
    """
    d = (coord.ra - cone_center.ra) ** 2 + (coord.dec - cone_center.dec) ** 2
    # The 0.0001 so we don't get edge effects
    return d < (cone_radius ** 2)

def find_ref(target: SkyCoord, refcoord: SkyCoord, threshold: u.degree = u.Quantity(0.05, u.deg)):
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


def photometry(x,y,aperture,sci_img):

    '''
    This function calculates the magnitude of any designated sources
    by conducting photometry on a masked and background subtracted image
    '''

    data = sci_img['SCI'].data
    gain = sci_img['SCI'].header['GAIN']
    exp_time = sci_img['SCI'].header['EXPTIME']

    # masking
    mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
    data = data.byteswap().newbyteorder()
    bkg = sep.Background(data, mask=mask)
    #imgdata_bkgsub = data - bkg.back()
    imgdata_bkgsub = data
    # create an error map (extra steps are taken to handle negative values)
    errmap = np.sqrt(np.sqrt(imgdata_bkgsub**2))/gain

    # Photometry
    source_pos = np.transpose((x, y))
    source_ap = CircularAperture(source_pos, r=aperture)
    source_flux = aperture_photometry(imgdata_bkgsub, source_ap, error=errmap)
    mag = -2.5 * np.log10(source_flux['aperture_sum'] * gain / exp_time)
    magerr = np.sqrt(((-5 / ((2 * np.log(10)) * (source_flux['aperture_sum'] * gain))) * (source_flux['aperture_sum_err'] * gain)) ** 2)

    return mag, magerr

# Imitate the layout of mags.py, but only doing instrumental photometry
# Using the ref catalog as a star catalog

files = glob('/home/pkotta/sdi_output_test/*.fits',recursive=True)
ims = [fits.open(f) for f in files]

#Science images:
sci_f = glob('/home/pkotta/GTAnd_SCI/*.fits',recursive=True)
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
ref_source = hdultocluster(ims, name="REF", tablename="ROBJ")
cat_source = hdultocluster(ims, name = 'CAT', tablename= 'OBJ')
var_source = hdultocluster(ims, name = 'XRT', tablename= 'XOBJ')
target_coords = SkyCoord(ims[0]['CAT'].data['ra'],ims[0]['CAT'].data['dec'], frame= 'icrs',unit='degree')

# Setup empty lists for storing instrumental magnitudes
inst_mag_li = []
inst_magerr_li = []

# TODO: placeholder for loop, don't want to go through and change indentation
for i in range(1):
    #comp_coords = find_ref(target_coord, refcoord) #in all images
    comp_coords = [refcoord[i] for i in np.random.randint(len(refcoord), size=15)]

    # Manually specifying stars of interest for now:
    #comp_coords = SkyCoord([11.25644924, 11.30799598, 11.32681067, 11.34687563, 11.35917106, 11.28148362], [41.59158802, 41.59097293, 41.59066171, 41.5905362, 41.59111281, 41.58777675], frame='icrs', unit='degree')
    #%%
    #Find the comp_stars in the ref hdu
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
    #next remove those indices in all the comp stars and sources, and remove those images
    comp_stars = [comp_stars[i][t_idx] for i in range(0,len(comp_stars))]
    comp_sources = [comp_sources[i][t_idx] for i in range(0,len(comp_sources))]
    ims_new = [ims[i] for i in t_idx[0]]
    #Now check for zero x, y, and a in the comp stars. If there are zero values, remove the comp source entirely
    c_idx = []
    for i in range(0,len(comp_sources)):
        if comp_sources[i]['x'].all()!=0 and comp_sources[i]['y'].all()!=0 and comp_sources[i]['a'].all()!=0:
            c_idx.append(i)
        else:
            print(comp_sources[i]['x'])
            print(comp_sources[i]['y'])
            print(comp_sources[i]['a'])
    comp_sources = np.array(comp_sources)[c_idx]
    comp_stars = np.array(comp_stars)[c_idx]

    for im in range(0,len(ims_new)):
        instrumental_mag = []
        in_magerr = []
        try:
            t = sci_ims[im][0].header['DATE-OBS']
            times = Time(t)
            #ts.append(times.mjd)
        except IndexError:
            t = sci_ims[im]['SCI'].header['OBSMJD']
            #ts.append(np.array(t))
        for i in range(0,len(comp_sources)):
            arr = photometry(comp_sources[i]['x'][im],comp_sources[i]['y'][im], comp_sources[i]['a'][im], sci_ims[im])
            instrumental_mag.append(arr[0][0])
            in_magerr.append(arr[1][0])
        inst_mag_li.append(instrumental_mag)
        inst_magerr_li.append(in_magerr)
    print(inst_mag_li)
    print(len(inst_mag_li[5]))
    print(len(inst_mag_li[25]))
    

# plot the instrumental magnitudes for each source across all avaliable images
fig, ax = plt.subplots(figsize=(10,7))
for i in range(len(inst_mag_li[0])):
    #ax.errorbar([l for l in range(len(inst_mag_li[0]))], [j[i] for j in inst_mag_li], yerr=[k[i] for k in inst_magerr_li])
    ax.plot([j[i] for j in inst_mag_li])
    #ax.errorbar([j[i] for j in inst_mag_li], yerr=[k[i] for k in inst_magerr_li])

plt.show()
