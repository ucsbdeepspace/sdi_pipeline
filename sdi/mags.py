"""
Photometry information generator

Programmed by Edgar Mao - 14 July 2022
"""

import click
from . import _cli as cli
import numpy as np
import matplotlib.pyplot as plt
import sep
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils.aperture import EllipticalAperture
from photutils.aperture import aperture_photometry
from glob import glob
from photutils.segmentation import make_source_mask
from astropy.time import Time
#Inter script functions
from .collate import hdultocluster, cluster_search_radec

# Cone Search Function
def _in_cone(coord: SkyCoord, cone_center: SkyCoord, cone_radius: u.degree):
    """
    Checks if SkyCoord coord is in the cone described by conecenter and
    cone_radius
    """
    d = (coord.ra - cone_center.ra) ** 2 + (coord.dec - cone_center.dec) ** 2
    # The 0.0001 so we don't get edge effects
    return d < (cone_radius ** 2)

# Photometry Function
def photometry(x,y,a,b,sci_img):
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
    imgdata_bkgsub = data - bkg.back()
    # create an error map (extra steps are taken to handle negative values)
    errmap = np.sqrt(np.sqrt(imgdata_bkgsub**2))/gain

    # Photometry
    source_pos = np.transpose((x, y))
    source_ap = EllipticalAperture(source_pos, a, b)
    source_flux = aperture_photometry(imgdata_bkgsub, source_ap, error=errmap)
    mag = -2.5 * np.log10(source_flux['aperture_sum'] * gain / exp_time)
    magerr = np.sqrt(((-5 / ((2 * np.log(10)) * (source_flux['aperture_sum'] * gain))) * (source_flux['aperture_sum_err'] * gain)) ** 2)

    return mag, magerr

# Reference Star Selection Function
def find_ref(target: SkyCoord, refcoord: SkyCoord, threshold: u.degree = u.Quantity(0.03, u.deg)):
    ref_comp = []
    # Select comparison stars around some target but exclude the target itself
    for coord in ref_coord:
        if _in_cone(target, coord, threshold) and not _in_cone(target, coord, u.Quantity(0.001, u.deg)):
            ref_comp.append(coord)
        else:
            pass
    return ref_comp

#Convert magnitudes of sources to SDSS g'
def norm(ref_table, cat_table):
    #gaia reference data
    ref_mag = ref_table['phot_g_mean_mag']
    g_rp_g = ref_table["phot_rp_mean_mag"]
    g_bp_g = ref_table["phot_bp_mean_mag"]

    gaia_mag_err = 0.16497
    r_err=0.066739
    #Making sure that all three values are present for calculating transformed mag
    gaia_mag_transformed = []
    r_transformed = []
    for i,j,k in zip(range(0,len(ref_mag)),range(0,len(g_rp_g)),range(0,len(g_bp_g))):
        if 25>ref_mag[i]>0.5  and 25>g_rp_g[j]>0.5  and 25>g_bp_g[k]>0.5:
            #transforming reference magnitudes to SDSS12 g
            gaia_mag_transformed.append(ref_mag[i]-0.13518+0.4625*(g_bp_g[k]-g_rp_g[j])+0.25171*(g_bp_g[k]-g_rp_g[j])**2-0.021349*(g_bp_g[k]-g_rp_g[j])**3)
            r_transformed.append(ref_mag[i]+0.12879-0.24662*(g_bp_g[k]-g_rp_g[j])+0.027464*(g_bp_g[k]-g_rp_g[j])**2+0.049465*(g_bp_g[k]-g_rp_g[j])**3)
        else:
            gaia_mag_transformed.append(np.nan)
            r_transformed.append(nan)

    gaia_mag_transformed = np.array(gaia_mag_transformed)
    r_transformed=np.array(r_transformed)

    return gaia_mag_transformed, r_transformed

def mags(hduls, read_ext="SCI", read_cat="XRT"):
    #TODO: right now these are not compatible with a single source
    #TODO: Need to make this entire function compatible with SUB hdu
    #Start by retrieveing radec for all reference stars
    #files = glob('/home/pkotta/sdi_output_test/*.fits',recursive=True)

    target_coords = []

    for hdul in hduls:
        cat = hdul[read_cat].data

        #Science images
        sci_ims = hdul[read_ext].data

        ref_ra = hdul['REF'].data['ra']
        ref_dec = hdul[0]['REF'].data['dec']
        refcoord = SkyCoord(ref_ra,ref_dec,frame = 'icrs',unit='degree')

        #Then use _in_cone to see which stars are close to the target star
        # Ex: For RR Lyrea:
        #target_coord = SkyCoord(11.291,41.508, frame = 'icrs', unit = 'degree')
        #comp_coords = find_ref(target_coord, refcoord) #in all images

        target_coords.append(SkyCoord(cat.data['ra'], cat.data['dec'], frame="icrs", unit='degree'))

    #Find the comp_stars in the ref hdu
    ref_source = hdultocluster(hduls, name="REF", tablename="ROBJ")
    cat_source = hdultocluster(hduls, name = 'CAT', tablename= 'OBJ')
    var_source = hdultocluster(hduls, name = 'XRT', tablename= 'XOBJ')
    for target_coord in target_coords:
        comp_coords = find_ref(target_coord, refcoord) #in all images
        var_sources = [cluster_search_radec(var_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords]
        var_coords = [SkyCoord(v['ra'],v['dec'], frame = 'icrs', unit = 'degree') for v in var_sources][0]
        # Check that comp_coords are not variable source from XRT
        nonvar_idx = []
        for v in var_coords:
            for c in range(0,len(comp_coords)):
                if v!=comp_coords[c]:
                    nonvar_idx.append(c)
            else:
                pass
        comp_coords_new = np.array(comp_coords)[list(set(nonvar_idx))]
        # Collate the comparision stars
        comp_sources = [cluster_search_radec(cat_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords_new]
        comp_stars = [cluster_search_radec(ref_source, coord.ra.deg, coord.dec.deg) for coord in comp_coords_new]
        #This is the target star in every image
        target = cluster_search_radec(cat_source, target_coord.ra.deg,target_coord.dec.deg)
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

        comp_sources = np.array(comp_sources)[c_idx]
        comp_stars = np.array(comp_stars)[c_idx]
        #For multiple reference stars
        mag_temp, r_temp = [norm(i,j)[0] for i,j in zip(comp_stars,comp_sources)]
        target_mag = []
        target_magerr = []
        ts = []
        for im in range(0, len(ims_new)):
            try:
                t = sci_ims[im][0].header['DATE-OBS']
                times = Time(t)
                ts.append(times.mjd)
            except IndexError:
                t = sci_ims[im]['SCI'].header['OBSMJD']
                ts.append(np.array(t))
            #For multiple ref stars
            #!TODO: Currently, comp_sources still have zeros? np.mean will not work anyway
            target_inst_mag = photometry(target['x'][im],target['y'][im], target['a'][im], sci_ims[im])[0][0]
            instrumental_mag = []
            in_magerr = []
            for i in range(0,len(comp_sources)):
                arr = photometry(comp_sources[i]['x'][im],comp_sources[i]['y'][im], comp_sources[i]['a'][im], sci_ims[im])
                instrumental_mag.append(arr[0][0])
                in_magerr.append(arr[1][0])

            bright_idx = []
            ref_magerr = 0.16497
            for i in range(0,len(instrumental_mag)):
                if instrumental_mag[i] < target_inst_mag: #20 is the lowest magnitude our pipeline can detect
                    bright_idx.append(i)
                else:
                    pass
            instrumental_mag = np.array(instrumental_mag)[bright_idx]
            ref_mag = np.array(mag_temp)[bright_idx]
            r_ref_mag = np.array(r_mag)[bright_idx]
            #Using 'a' as a sloppy alternative to aperture. Maybe look into sep.kron_radius or flux_radius. aperture is a radius.
            #Have to select index [0] for each mag to get just the numerical value without the description of the column object


            #Performing the linear fit
            #First find all places where there are non-nan values:
            idx = np.where([np.isnan(r)==False for r in ref_mag])[0]
            x = [instrumental_mag[i] for i in idx]
            y = [ref_mag[i] for i in idx]
            r = [r_ref_mag[i] for i in idx]
            
            def fit_eqn(flux, zp, a1, a2, a3,  a4, a5, t, t_m, exp_time):#flux is the x axis
                -2.5*np.log10(flux)+zp = -g+a1*(g-r)+a2*air_mass+a3*air_mass*(g-r)+a4*(t-t_m)+a5*(t-t_m)**2 - 2.5*np.log10(exp_time)
            
            from scipy.optimize import curve_fit
            
            fit, sum_sq_resid, rank, singular_values, rcond = np.polyfit(x[1:], y[1:], 1, full=True)
            fit_fn = np.poly1d(fit)
            residuals = fit_fn(x)-y
            print(len(y))
            '''
            fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
            plt.plot(x,y, 'yo', x, fit_fn(x), '--k')
            plt.xlabel('instrumental mag')
            plt.ylabel('reference magnitude')
            plt.show()
            '''
            # fit_fn is now a function which takes in x and returns an estimate for y,
            #Use the fit from above to calculate the target magnitude
            #target_inst_mag = photometry(target['x'][im],target['y'][im], target['a'][im], sci_ims[im])[0][0]
            target_mag.append(fit_fn(target_inst_mag))
            target_magerr.append(np.std(residuals))
        #error propagation
        error = []
        for i in range(len(target_magerr)):
            if np.abs(target_magerr[i])>np.mean(target_magerr)+np.std(target_magerr):
                error.append(target_magerr[i])
            else:
                error.append(in_magerr[i])

        rows = zip(target_mag,target_magerr, ts)
        print('Done!')
        '''
        import csv
        wfname = 'var_bright_test'+str(len(str(target_coord)))+'.csv'
        with open(wfname, "w") as f:
            writer = csv.writer(f)
            for row in rows:
                writer.writerow(row)
        '''
        new_hdu = fits.BinTableHDU(rows, header=None, name="MAG")
        hdul.append(new_hdu)
        yield hdul


@cli.cli.command("mags")
@click.option("-i", "--image", default="SCI", help="The image to do photometry on.")
@click.option("-c", "--catalog", default="XRT", help="The catalog that the function is referencing the source to.")
@cli.operator

## mags function wrapper

def mags_cmd(hduls, image="SCI", catalog="XRT"):
    """
    Generates photometry information for designated sources.
    """
    return mags(hduls, image, catalog)
