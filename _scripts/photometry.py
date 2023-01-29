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
from astropy.time import Time
from collate import hdultocluster, cluster_search_radec
import csv
import warnings
warnings.filterwarnings(action='ignore')


def _in_cone(coord: SkyCoord, cone_center: SkyCoord, cone_radius: u.degree):
    """
    Checks if SkyCoord coord is in the cone described by conecenter and
    cone_radius
    """
    d = (coord.ra - cone_center.ra) ** 2 + (coord.dec - cone_center.dec) ** 2
    return d < (cone_radius ** 2)


def find_ref(target: SkyCoord, ref_coords: SkyCoord, threshold: u.degree = u.Quantity(0.06, u.deg)):
    ref_comp = []
    for coord in ref_coords:
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

def norm(ref_table):
    ref_mag = ref_table['phot_g_mean_mag']
    g_rp_g = ref_table["phot_rp_mean_mag"]
    g_bp_g = ref_table["phot_bp_mean_mag"]

    #Making sure that all three values are present for calculating transformed mag
    gaia_mag_transformed = []
    r_transformed = []
    diffs = []
    for i,j,k in zip(range(0,len(ref_mag)),range(0,len(g_rp_g)),range(0,len(g_bp_g))):
        if 25 > ref_mag[i] > 0.5  and 25 > g_rp_g[j] > 0.5  and 25 > g_bp_g[k] > 0.5:
            #transforming reference magnitudes to SDSS12 g
            diffs.append(g_bp_g[k]-g_rp_g[j])
            gaia_mag_transformed.append(ref_mag[i] - 0.13518 + 0.4625*(g_bp_g[k] - g_rp_g[j]) + 0.25171 * (g_bp_g[k] - g_rp_g[j])**2-0.021349 * (g_bp_g[k]-g_rp_g[j])**3)
            r_transformed.append(ref_mag[i] + 0.12879 - 0.24662 * (g_bp_g[k] - g_rp_g[j]) + 0.027464 * (g_bp_g[k]-g_rp_g[j])**2 + 0.049465*(g_bp_g[k]-g_rp_g[j])**3)
        else:
            gaia_mag_transformed.append(np.nan)
            r_transformed.append(np.nan)
    gaia_mag_transformed = np.array(gaia_mag_transformed)
    r_transformed=np.array(r_transformed)

    return gaia_mag_transformed, r_transformed

def photometry(x, y, aperture, im):

    '''
    This function calculates the magnitude of any designated sources
    by conducting photometry on a masked and background subtracted image
    '''

    data = im['SCI'].data
    gain = im['SCI'].header['GAIN']

    # masking
    mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
    data = data.byteswap().newbyteorder()
    bkg = sep.Background(data, mask=mask)
    imgdata_bkgsub = data - bkg.back()
    errmap = np.sqrt(np.sqrt(imgdata_bkgsub**2))/gain

    # Photometry
    source_pos = np.transpose((x, y))
    source_ap = CircularAperture(source_pos, r=aperture)
    source_flux = aperture_photometry(imgdata_bkgsub, source_ap, error=errmap)
    mag = -2.5 * np.log10(source_flux['aperture_sum'] * gain)
    magerr = np.sqrt(((-5 / ((2 * np.log10(10)) * (source_flux['aperture_sum'] * gain))) * (source_flux['aperture_sum_err'] * gain)) ** 2)

    return mag, magerr


def Error_Finder(bkgd, mag, gain, FWHM):
    cnst_err = gain * bkgd * np.pi * (FWHM/2)**2
    mag_err = 1.086/gain**2 * np.sqrt(cnst_err + gain*10**(-0.4*mag))/10**(-0.4*mag)
    return mag_err


def bkgd_fit(bkgd_guess, zp, i_mag, i_g, gain, fwhm): #Recursively, minimizes the median absolute deviation of the function
    params = []
    median_list = []
    for bkgd in bkgd_guess:
        params.append(bkgd) #append this guess set to the list.
        model_output = Error_Finder(bkgd, i_mag, gain, fwhm)
        abs_devs = np.abs(i_g - model_output + zp)
        med_abs_dev = np.median(abs_devs)
        median_list.append(med_abs_dev)
    best_params = params[np.where(median_list == np.min(median_list))[0][0]]
    return best_params


def mad_curve_fit(ZP_guess, c0_guesses, i_mag, g, g_r): #Recursively, minimizes the median absolute deviation of the function
    params = []
    median_list = []
    for zp in ZP_guess:
            for c0 in c0_guesses:
                guess_params = [zp, c0] #make an array of guesses
                params.append(guess_params) #append this guess set to the list.
                abs_devs = np.abs(g - i_mag - zp - c0 * g_r)
                med_abs_dev = np.median(abs_devs)
                median_list.append(med_abs_dev)
    best_params = params[np.where(median_list == np.min(median_list))[0][0]]
    return best_params

def get_reference_stars(cat_catalog, ref_catalog, variable_catalog, ref_coords):
    #Remove reference stars that are also variable
    variable_stars_in_catalog = [cluster_search_radec(variable_catalog, coord.ra.deg, coord.dec.deg) for coord in ref_coords]
    variable_star_coords = [SkyCoord(v['ra'],v['dec'], frame = 'icrs', unit = 'degree') for v in variable_stars_in_catalog][0]
    nonvar_idx = []
    for v in variable_star_coords:
        for c in range(0, len(ref_coords)):
            if v!=ref_coords[c]:
                nonvar_idx.append(c)
            else:
                pass
    ref_coords = np.array(ref_coords)[list(set(nonvar_idx))]
    
    #These are the comp stars in the cat table itself. so these stars are in our images
    sdi_reference_stars = [cluster_search_radec(cat_catalog, coord.ra.deg, coord.dec.deg) for coord in ref_coords]
    #And these are the same refence stars but from gaia
    ref_stars = [cluster_search_radec(ref_catalog, coord.ra.deg, coord.dec.deg) for coord in ref_coords]
    
    #Now check for zero x, y, and a in the comp stars. If there are zero values, remove the comp source entirely
    c_idx = []
    for i in range(0,len(sdi_reference_stars)):
        if sdi_reference_stars[i]['x'].all() !=0 and sdi_reference_stars[i]['y'].all() !=0 and sdi_reference_stars[i]['a'].all()!=0:
            c_idx.append(i)
    
    sdi_reference_stars = np.array(sdi_reference_stars)[c_idx] #The flux of the reference stars in our images
    ref_stars = np.array(ref_stars)[c_idx] #The apparent magnitude of reference stars in Gaia
    return sdi_reference_stars, ref_stars


#-----------Retrieve files---------------
files = glob(r'C:\Users\Sam Whitebook\Documents\Visual Studio 2010\Projects\Lubin Lab\Light_Curves\sdi_output\*.fits', recursive=True)
ims = [fits.open(f) for f in files]

#Science images:
sci_f = glob(r'C:\Users\Sam Whitebook\Documents\Visual Studio 2010\Projects\Lubin Lab\Light_Curves\GTAnd_SCI\*.fits', recursive=True)
sci_ims = [fits.open(f) for f in sci_f]

del files, sci_f


#----------Collect our sources-----------------
ref_ra = ims[0]['REF'].data['ra']
ref_dec = ims[0]['REF'].data['dec']
ref_coords = SkyCoord(ref_ra,ref_dec,frame = 'icrs',unit='degree')

ref_catalog = hdultocluster(ims, name="REF", tablename="ROBJ") #Reference stars from Gaia DR2
cat_catalog = hdultocluster(ims, name = 'CAT', tablename= 'OBJ') #All sources found in the image
variable_catalog = hdultocluster(ims, name = 'XRT', tablename= 'XOBJ') #All variable sources found after subtraction

target_coord = SkyCoord(11.291, 41.508, frame = 'icrs', unit = 'degree') 

target_star_in_catalog = cluster_search_radec(cat_catalog, target_coord.ra.deg, target_coord.dec.deg)

del ref_ra, ref_dec


#----------------Find reference stars close to the target------------------
sdi_reference_stars, ref_stars = get_reference_stars(cat_catalog, ref_catalog, variable_catalog, ref_coords)


#This is the target star in every image
target_star = cluster_search_radec(cat_catalog, target_coord.ra.deg, target_coord.dec.deg)

#Remove all the images for which the target star could not be found by removing places that have zeroes
target_not_present_idx = np.where(target_star['x'] != 0)
target = target_star[target_not_present_idx]

#next remove those indices in all the comp stars and sources, and remove those images
ref_stars = [ref_stars[i][target_not_present_idx] for i in range(0,len(ref_stars))]
sdi_reference_stars = [sdi_reference_stars[i][target_not_present_idx] for i in range(0,len(sdi_reference_stars))]
ims_with_target = [i for i in target_not_present_idx[0]]

#Convert Gaia g, r magnitudes to SDSS g' r' magnitudes. This is because our images are taken in SDSS g' r' filters
ref_stars = np.array([norm(i) for i in ref_stars])
reference_star_g_mag = ref_stars[:,0,0]
reference_star_r_mag = ref_stars[:,1,0]

del target_coord, ref_coords, ref_catalog, cat_catalog, variable_catalog, ims

#----------------Estimate Uncertainty and Target Magnitude------------------
target_mag = []
target_magerr = []
times = []
for im in ims_with_target:
    try:
        t = sci_ims[im][0].header['DATE-OBS']
        time = Time(t)
        t = time.mjd
        times=times.append(t)
    except IndexError:
        t = sci_ims[im]['SCI'].header['OBSMJD']
        times.append(np.array(t))
    
    gain = sci_ims[im][0].header['GAIN']
  
    #Calculate the instrumental magnitude (flux) of the star
    #Using 'a' as a sloppy alternative to aperture. Maybe look into sep.kron_radius or flux_radius. aperture is a radius.
    #Have to select index [0] for each mag to get just the numerical value without the description of the column object
    target_inst_mag = photometry(target['x'][im],target['y'][im], target['a'][im], sci_ims[im])[0][0]

    instrumental_mag = []
    in_magerr = []
    for i in range(0,len(sdi_reference_stars)):
        arr = photometry(sdi_reference_stars[i]['x'][im],sdi_reference_stars[i]['y'][im], sdi_reference_stars[i]['a'][im], sci_ims[im])
        instrumental_mag.append(arr[0][0])
        in_magerr.append(arr[1][0])

    #Only select reference stars that are brighter than 20th magnitude
    bright_stars_idx = []
    ref_magerr = 0.16497
    for i in range(0,len(instrumental_mag)):
        if reference_star_g_mag[i] < 18: #20 is the lowest magnitude our pipeline can detect
            bright_stars_idx.append(i)
        else:
            pass
    instrumental_mag = np.array(instrumental_mag)[bright_stars_idx]
    ref_mag = np.array(reference_star_g_mag)[bright_stars_idx]
    r_ref_mag = np.array(reference_star_r_mag)[bright_stars_idx]
    
    #First find all places where there are non-nan values:
    nan_idx = np.where([np.isnan(r)==False for r in ref_mag])[0]
    instrumental_mag = [instrumental_mag[i] for i in nan_idx]
    g = [ref_mag[i] for i in nan_idx]
    r = [r_ref_mag[i] for i in nan_idx]
    instrumental_mag = np.array(instrumental_mag)
    g = np.array(g)
    r = np.array(r)
    g_r = g - r
    inst_minus_g = instrumental_mag - g

    #Assuming Gaia's magnitudes are the true values of the star, we can use our instrumental mags in g to estimate uncertainty.
    zp_guesses = np.linspace(25, 30, 80)
    bkgd_guesses = np.linspace(15000, 25000, 100)
    c0_guesses = np.linspace(-1.5, 1.5, 20)

    best_params = mad_curve_fit(zp_guesses, c0_guesses, instrumental_mag, g, g_r)
    zp = best_params[0]
    c0 = best_params[1]

    bkgd = bkgd_fit(bkgd_guesses, zp, instrumental_mag, inst_minus_g, 1.6, 12.5)
    
    target_mag_err = Error_Finder(bkgd, target_inst_mag, 1.6, 12.5)

    target_mag.append(target_inst_mag + zp)
    target_magerr.append(target_mag_err)

rows = zip(target_mag, target_magerr, times)
wfname = 'GTAnd_least_abs_dev_g_i.csv'
with open(wfname, 'w') as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)