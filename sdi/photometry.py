import glob
import sep
import astropy.units as u
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import numpy as np
import astroalign as aa
from astropy.io import fits
from astroquery.ipac.irsa import Irsa
from astroquery.sdss import SDSS
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy.coordinates import SkyCoord
from matplotlib.colors import LogNorm
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus, ApertureStats
from regions import CirclePixelRegion, PixCoord
from photutils import centroids
import multiprocessing as mp
from multiprocessing import shared_memory
from . import _cli as cli
import click

def photometry(hduls, name, directory):
    os.chdir(directory)
    data = [h[name].data for h in hduls]
    template = np.median(data, axis = 0)
    bkg_phot = sep.Background(template)
    extracted_phot = sep.extract(template - bkg_phot.back())
    w = WCS(hdus[0][1].header)
    print(len(extracted_phot))
    print(w)
    sources = []
    for c, src in enumerate(extracted_phot): #indexes extracted sources by try number to find reference stars
        x = src['x']
        y = src['y']
        coord = pixel_to_skycoord(x, y, w).transform_to('icrs') #gives wcs transformation for pixel coordinates
        search = SDSS.query_crossid(coord, fields = ['ra', 'dec', 'psfMag_r', 'psfMagErr_r'], radius = 15 * u.arcsec, region = False) #narrow field cone search to find source based on ra, dec.
        if search: #if SDSS query returned results, continue
            if search['psfMag_r'] < 16 and search['type'] == 'STAR': #filters search results by sources that are brighter than magnitude 16, and of type star.
                ref = dict(ra = search['ra'], dec = search['dec'], x_coord = x, y_coord = y, x_min = src['xmin'], x_max = src['xmax'],rad = None, g_mag = search['psfMag_r'], g_mag_err = search['psfMagErr_r'], source_id = c, calibrated_mags = [], instrumental_mags = [], inst_mag_errs= []) #creates dictionary item with source position (ra, dec, x, y), source extent, and mags and errors as reported by SDSS.
                sources.append(ref)
            else:
                ref = dict(ra = coord, dec = coord, x_coord = x, y_coord = y, x_min = src['xmin'], x_max = src['xmax'],rad = None, g_mag = None, g_mag_err = None, source_id = c, calibrated_mags = [], instrumental_mags = [], inst_mag_errs = []) #if the search returns results but doesn't meet criteria, treat same as if no search was returned. Cleaner way to do this?
                sources.append(ref)
        else:
            ref = dict(ra = coord, dec = coord, x_coord = x, y_coord = y, x_min = src['xmin'], x_max = src['xmax'], rad = None, g_mag = None, g_mag_err = None, source_id = c, calibrated_mags = [],instrumental_mags = [], inst_mag_errs = []) #sources that are not reference stars simply will not have SDSS magnitudes.
            sources.append(ref)
    for i, image in enumerate(aligned):
        print(i)
        N_r = hdus[i][1].header["RDNOISE"] #readout noise
        for source in sources:
            coords = [source['x_coord'], source['y_coord']]
            pcoords = PixCoord(source['x_coord'], source['y_coord'])
            radius_i = (source['x_max'] - source['x_min']) / 2
            source['rad'] = (radius_i)
            radius_o_0 = radius_i + 5
            radius_o_1 = radius_o_0 + 5

            source_circle = CirclePixelRegion(pcoords, radius_i).to_mask() #makes region of source shape
            source_aperture = source_circle.cutout(image) #gets data of source

            background_annulus = CircularAnnulus(coords, radius_o_0, radius_o_1)
            background_mean = ApertureStats(image, background_annulus).mean

            source_flux_pix = source_aperture-(source_circle*background_mean) #pixel wise background subtraction
            source_flux_total = np.sum(source_flux_pix)

            readout_sum_square = np.sum(source_circle*np.float64(N_r**2)) #applies square readout noise to source array shape, then adds. Gives sum of square readout noise over back subtracted source.

            delta_n = (readout_sum_square + source_flux_total + (((radius_i**2)/((radius_o_1**2)-(radius_o_0**2)))**2)*(readout_sum_square + aperture_photometry(image, background_annulus)['aperture_sum'][0]))**(1/2) #this is the stuff for SNR

            if source_flux_total <= 0:
                inst_mag = -2.5*np.log10(abs(source_flux_total)) # For now, the case where the background is oversubtracted from LCO is handled in this way but this is probably not the correct way to do this.
                delta_m = 2.5*np.log10(np.e)*abs(delta_n/source_flux_total)
            else:
                inst_mag = -2.5*np.log10(source_flux_total)
                delta_m = 2.5*np.log10(np.e)*abs(delta_n/source_flux_total)
            source['instrumental_mags'].append(inst_mag)
            source['inst_mag_errs'].append(delta_m)

    #building calibration model:
    res = []
    inst_mags = [np.mean(source['instrumental_mags']) for source in sources if source['g_mag'] != None]
    sky_mags = [source['g_mag'] for source in sources if source['g_mag'] != None]
    print(len(inst_mags), len(sky_mags))

    #Makes linear model for calibration:
    #This is the first round of modeling, with outliers.
    p0 = np.polyfit(inst_mags, sky_mags, deg = 1)
    x = np.arange(-15, 0)
    y = p0[0]*x + p0[1]
    plt.plot(x, y, color = 'b', label = "Model With Outliers")
    diffs = [s['g_mag']- (np.mean(s['instrumental_mags'])*p0[0] + p0[1]) for s in sources if s['g_mag'] != None]
    stdv = np.std(diffs)

    inst_mags_final = []
    sky_mags_final = []
    outlier_inst =[]
    outlier_sky =[]

    for diff in diffs: #rudementary sigma clipping to remove outliers from calibration model.
        if diff < stdv:
            i = diffs.index(diff)
            inst_mags_final.append(inst_mags[i])
            sky_mags_final.append(sky_mags[i])
        else:
            i = diffs.index(diff)
            outlier_inst.append(inst_mags[i])
            outlier_sky.append(sky_mags[i])
    p1 = np.polyfit(inst_mags_final, sky_mags_final, deg = 1) #recalculates calibration model without outliers.
    #p2 = np.polyfit(inst_mags_final, sky_mags_final, deg = 0)
    #print(p2[0])
    print("first try: {}".format(p0)) #prints slopes of each model. In theory, they should come out to around 1.
    print("second try: {}".format(p1))

    plt.scatter(outlier_inst, outlier_sky, color = 'b', label = "Outliers")
    plt.scatter(inst_mags_final, sky_mags_final, color = 'r', label = "Kept")
    plt.plot(x, [i*p1[0] + p1[1] for i in x], color = 'r', label = "Model Without Outliers")
    plt.plot(x, [i+ p1[1] for i in x], color = 'g', label = "unity")
    plt.xlabel("Instrumental Magnitude SDSS g-band")
    plt.ylabel("SDSS Reference Magnitude g-band")
    plt.title("Instrumental vs Reference Magnitude")
    plt.legend()
    #plt.savefig("F:/SDI/Section32Figures/calibrationplot.png", dpi = 1000)
    #plt.savefig("/Users/lucaangeleri/Documents/LCO/sec17figures/calibrationplot.png", dpi = 1000)
    #plt.show()
    plt.savefig('calibration with reference stars.png')

    #add calibrated mags to sources:
    for source in sources:
        vals = []
        for val in source['instrumental_mags']:
            cal = np.float64(val*p1[0] + p1[1])
            vals.append(cal)
        source['calibrated_mags'] = vals #probably a cleaner way to do this part but was having issue where calibrated magnitudes were being added to dict as individual arrays

    for source in sources:
        r =  np.arange(0, len(source['calibrated_mags']), 1)
        plt.errorbar(r, source['calibrated_mags'], yerr = source['inst_mag_errs'], linestyle = 'none', marker = 'o', color = 'b')

        Chis = []
        avg_mag = np.mean(source['calibrated_mags'])
        for i, m in enumerate(source['calibrated_mags']):
            chi_i = ((m- avg_mag)**2)/(source['inst_mag_errs'][i]**2)
            Chis.append(chi_i)
        dof = len(source['calibrated_mags']) - 1
        chi_dof = np.sum(Chis)/dof
        #plt.title("X, Y: {}, {}; RA, DEC: {}, {}, CHI2: {}, ID: {}".format("%.2f" % source['x_coord'],"%.2f" % source['y_coord'],"%.4f" % source['ra'][0],"%.4f" % source['dec'][0], "%.2f" %chi_dof, source['source_id'] )) #get only to display a few decimal places
        plt.title("X, Y: {}, {}; CHI2: {}, ID: {}".format("%.2f" % source['x_coord'],"%.2f" % source['y_coord'], "%.2f" %chi_dof, source['source_id'] ))
        plt.plot(r, np.ones(len(r))*avg_mag, label = "TRIPP AVG MAG", linestyle = '--', color = 'b')
        if source['g_mag'] != None:
            plt.plot(r, np.ones(len(r))*source['g_mag'], linestyle = '--', color = 'r', label = "SDSS REPORTED MAG" )
        plt.legend()
        plt.savefig('Source {} Lightcurve.png'.format(source['source_id'], dpi = 500))

        plt.title("Source Number: {}, Position: {}, {}".format(source['source_id'], source['x_coord'], source['y_coord']))
        circle0 = plt.Circle((source['x_coord'], source['y_coord']), source['rad'], color = 'r', linewidth= .25, fill = False)
        circle1 = plt.Circle((source['x_coord'], source['y_coord']), source['rad'] + 5, color = 'b', linewidth= .25,fill = False)
        circle2 = plt.Circle((source['x_coord'], source['y_coord']), source['rad'] +5, color = 'g', linewidth= .25,fill = False)
        ax = plt.gca()
        ax.add_patch(circle0)
        ax.add_patch(circle1)
        ax.add_patch(circle2)
        plt.imshow(template, cmap = 'gray', norm = LogNorm(vmin = 1, vmax = 200), origin='lower')
        plt.savefig("Source location {}.png".format(source['source_id']), format = 'png', dpi = 500)



@cli.cli.command("photometry")
@click.option('-d', '--directory', type=str, help="Specify path to directory to save fitsfiles.", default="./")
@click.option("-n", "--name", default="SUB", help="The HDU to run photometry on")
@cli.operator
## Photometry function wrapper
def photometry_cmd(hduls, name='SUB', directory = './'):
    """
    Returns lightcurves from a set of subtracted or unsubtracted images from a template image\n
    Arguments:\n
        hduls -- list of fits hduls\n
    """
    return photometry(hduls, name, directory)
