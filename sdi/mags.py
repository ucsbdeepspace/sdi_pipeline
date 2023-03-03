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
import warnings
import scipy as sp
import csv
from scipy import special
from scipy.optimize import *

def _in_cone(coord: SkyCoord, cone_center: SkyCoord, cone_radius: u.degree):
    """
    Checks if SkyCoord coord is in the cone described by conecenter and
    cone_radius
    """
    d = (coord.ra - cone_center.ra) ** 2 + (coord.dec - cone_center.dec) ** 2
    # The 0.0001 so we don't get edge effects
    return d < (cone_radius ** 2)

def find_ref(target: SkyCoord, refcoord: SkyCoord, threshold: u.degree = u.Quantity(0.06, u.deg)):
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
    #gaia reference data
    ref_mag = ref_table['phot_g_mean_mag']
    g_rp_g = ref_table["phot_rp_mean_mag"]
    g_bp_g = ref_table["phot_bp_mean_mag"]

    gaia_mag_err = 0.16497
    r_err=0.066739
    #Making sure that all three values are present for calculating transformed mag
    gaia_mag_transformed = []
    r_transformed = []
    diffs = []
    for i,j,k in zip(range(0,len(ref_mag)),range(0,len(g_rp_g)),range(0,len(g_bp_g))):
        if 25>ref_mag[i]>0.5  and 25>g_rp_g[j]>0.5  and 25>g_bp_g[k]>0.5:
            #transforming reference magnitudes to SDSS12 g
            diffs.append(g_bp_g[k]-g_rp_g[j])
            gaia_mag_transformed.append(ref_mag[i]-0.13518+0.4625*(g_bp_g[k]-g_rp_g[j])+0.25171*(g_bp_g[k]-g_rp_g[j])**2-0.021349*(g_bp_g[k]-g_rp_g[j])**3)
            r_transformed.append(ref_mag[i]+0.12879-0.24662*(g_bp_g[k]-g_rp_g[j])+0.027464*(g_bp_g[k]-g_rp_g[j])**2+0.049465*(g_bp_g[k]-g_rp_g[j])**3)
        else:
            gaia_mag_transformed.append(np.nan)
            r_transformed.append(np.nan)
    gaia_mag_transformed = np.array(gaia_mag_transformed)
    r_transformed=np.array(r_transformed)

    return gaia_mag_transformed, r_transformed

def photometry(x,y,aperture,sci_img):

    '''
    This function calculates the magnitude of any designated sources
    by conducting photometry on a masked and background subtracted image
    '''

    data = sci_img['SCI'].data
    gain = sci_img['SCI'].header['GAIN']
    #exp_time = sci_img['SCI'].header['EXPTIME']

    # masking
    mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
    data = data.byteswap().newbyteorder()
    bkg = sep.Background(data, mask=mask)
    imgdata_bkgsub = data - bkg.back()
    # create an error map (extra steps are taken to handle negative values)
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

def chi_sq(params, x, f, err):
    m = params[0]
    b = params[1]
    chi = np.sum((f - m*x - b)**2/err**2)
    return chi

def mad_fit(ZP_guess, bkgd_guess, i_mag, i_g, gain, fwhm): #Recursively, minimizes the median absolute deviation of the function
    params = []
    median_list = []
    for zp in ZP_guess:
        for bkgd in bkgd_guess:
                guess_params = [zp, bkgd] #make an array of guesses
                params.append(guess_params) #append this guess set to the list.
                model_output = Error_Finder(bkgd, i_mag, gain, fwhm)
                abs_devs = np.abs(i_g - g_r - model_output + zp)
                med_abs_dev = np.median(abs_devs)
                median_list.append(med_abs_dev)
    best_params = params[np.where(median_list == np.min(median_list))[0][0]]
    return best_params

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

def mags(hduls, read_ext):
    """
    Function for generating the apparent magnitude of any arbitrary sources
    """
    # Retrieve all images
    ims = [hdul for i in hduls]
    # Collect sources
    refcoord = SkyCoord(ref_ra,ref_dec,frame = 'icrs',unit='degree')

    #Eachof these tables are the one source across all images
    ref_source = hdultocluster(ims, name="REF", tablename="ROBJ") #Reference stars from Gaia DR2
    cat_source = hdultocluster(ims, name = 'CAT', tablename= 'OBJ') #All sources found in the image
    var_source = hdultocluster(ims, name = read_ext, tablename= 'XOBJ') #All variable sources found after subtraction
    target_coords = SkyCoord(ims[0]['CAT'].data['ra'],ims[0]['CAT'].data['dec'], frame= 'icrs',unit='degree')
# print(target_coords)
    #Test source RR Lyra
    #target_coord = SkyCoord(11.291,41.508, frame = 'icrs', unit = 'degree') #RR Lyrae
    # target_coord = SkyCoord(11.35917106, 41.59111281, frame = 'icrs', unit = 'degree') #Some other random star
    target = cluster_search_radec(cat_source, target_coords.ra.deg, target_coords.dec.deg)

    # Find reference stars
    compcoords = refcoord
    # Remove variable stars in reference stars catalog
    var_sources = [cluster_search_radec(var_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords]
    var_coords = [SkyCoord(v['ra'],v['dec'], frame = 'icrs', unit = 'degree') for v in var_sources][0]

    nonvar_idx = []
    for v in var_coords:
        for c in range(0,len(comp_coords)):
            if v!=comp_coords[c]:
                nonvar_idx.append(c)
            else:
                pass
    comp_coords_new = np.array(comp_coords)[list(set(nonvar_idx))]

    comp_sources = [cluster_search_radec(cat_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords_new]
    comp_stars = [cluster_search_radec(ref_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords_new]
    
    # Remove all the images for which the target star could not be found by removing places that have zeroes
    t_idx = np.where(target['x']!=0)
    target = target[t_idx]

    # next remove those indices in all the comp stars and sources, and remove those images
    comp_stars = [comp_stars[i][t_idx] for i in range(0,len(comp_stars))]
    comp_sources = [comp_sources[i][t_idx] for i in range(0,len(comp_sources))]
    ims_new = [i for i in t_idx[0]]

    # Now check for zero x, y, and a in the comp stars. If there are zero values, remove the comp source entirely
    c_idx = []
    for i in range(0,len(comp_sources)):
        if comp_sources[i]['x'].all()!=0 and comp_sources[i]['y'].all()!=0 and comp_sources[i]['a'].all()!=0:
            c_idx.append(i)

    comp_sources = np.array(comp_sources)[c_idx] # The flux of the reference stars in our images
    comp_stars = np.array(comp_stars)[c_idx] # The apparent magnitude of reference stars in Gaia

    # Convert Gaia g, r magnitudes to SDSS g' r' magnitudes. This is because our images are taken in SDSS g' r' filters
    norm_compstar = np.array([norm(i) for i in comp_stars])
    mag_temp = norm_compstar[:,0,0]
    r_temp = norm_compstar[:,1,0]
    target_mag = []
    target_magerr = []
    ts = []
    times = [im[0].header['OBSMJD'] for im in sci_ims]
    for im in ims_new:
        try:
            t = sci_ims[im][0].header['DATE-OBS']
            times = Time(t)
            time = times.mjd
            ts.append(time)
        except IndexError:
            t = sci_ims[im]['SCI'].header['OBSMJD']
            ts.append(np.array(t))

        gain = sci_ims[im][0].header['GAIN']
    #Calculate the instrumental magnitude (flux) of the star
        target_inst_mag = photometry(target['x'][im],target['y'][im], target['a'][im], sci_ims[im])[0][0]
        target_mag_err = photometry(target['x'][im],target['y'][im], target['a'][im], sci_ims[im])[1][0]
        instrumental_mag = []
        in_magerr = []
        for i in range(0,len(comp_sources)):
            arr = photometry(comp_sources[i]['x'][im],comp_sources[i]['y'][im], comp_sources[i]['a'][im], sci_ims[im])
            instrumental_mag.append(arr[0][0])
            in_magerr.append(arr[1][0])

        #Only select reference stars that are brighter than 20th magnitude
        bright_idx = []
        ref_magerr = 0.16497
        for i in range(0,len(instrumental_mag)):
            if mag_temp[i] < 18: #20 is the lowest magnitude our pipeline can detect
                bright_idx.append(i)
            else:
                pass
        instrumental_mag = np.array(instrumental_mag)[bright_idx]
        ref_mag = np.array(mag_temp)[bright_idx]
        r_ref_mag = np.array(r_temp)[bright_idx]
    
        #First find all places where there are non-nan values:
        idx = np.where([np.isnan(r)==False for r in ref_mag])[0]
        i_mag = [instrumental_mag[i] for i in idx]
        g = [ref_mag[i] for i in idx]
        r = [r_ref_mag[i] for i in idx]
        i_mag = np.array(i_mag)
        g = np.array(g)
        r = np.array(r)
        g_r = g-r
        i_g = i_mag - g

        #Assuming Gaia's magnitudes are the true values of the star, we can use our instrumental mags in g to estimate uncertainty.
        zp_guesses = np.linspace(25, 30, 80)
        bkgd_guesses = np.linspace(15000, 25000, 100)
        c0_guesses = np.linspace(-1.5, 1.5, 20)

        best_params = mad_curve_fit(zp_guesses, c0_guesses, i_mag, g, g_r)
        zp = best_params[0]
        c0 = best_params[1]

        bkgd = bkgd_fit(bkgd_guesses, zp, i_mag, i_g, 1.6, 12.5)
        
        #     #===============scipy g-r fit=====================
        # def fit(g_inst, zp, c0):
        #     g = g_inst + zp + c0*g_r
        #     return g
    
        target_mag_err = Error_Finder(bkgd, target_inst_mag, 1.6, 12.5)
        # gr_popt, gr_pcov = sp.optimize.curve_fit(fit, i_mag, g, p0 = [0, 23], sigma=mag_err)
        # _, zp_err = np.sqrt(np.diag(gr_pcov))
        # print('zp = ', gr_popt[0], 'c0= ', gr_popt[1])

        target_mag.append(target_inst_mag + zp)
        target_magerr.append(target_mag_err)

    # Output results (the output at the moment is csv files)
    rows = zip(target_mag, target_magerr, ts)
    wfname = 'source.csv'
    with open(wfname, 'w') as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)


@cli.cli.command("mags")
@click.option("-r", "--read_ext", default="XRT",
              help="An index number or ext name that identifies the target sources you want.") 
@cli.operator

def mags_cmd(hduls, threshold, read_ext, write_ext):
    """
    Uses sep to find sources in ImageHDU data.
    """
    try:
        read_ext = int(read_ext)
    except ValueError:
        pass
    return mags(hduls, read_ext)
