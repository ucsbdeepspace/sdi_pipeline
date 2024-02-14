# import csv
# import glob
import sep
import astroalign as aa
import astropy.units as u
import astropy.coordinates as coord
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.io import fits
from astroquery.ipac.irsa import Irsa
from astroquery.sdss import SDSS
import concurrent.futures
import click
from . import _cli as cli
from .combine import combine
import matplotlib.pyplot as plt
import numpy as np
from scipy import odr
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm, TABLEAU_COLORS, CSS4_COLORS
import multiprocessing as mp
from multiprocessing import shared_memory
import os
from regions import CirclePixelRegion, PixCoord
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus, ApertureStats
from photutils.aperture import CircularAperture, aperture_photometry
from photutils import background as bck
from photutils.segmentation import detect_sources,detect_threshold,SourceCatalog
from photutils import detection
from photutils import centroids
import time


class Night: #Each observing night is initialized as a Night object
    def __init__(self, img_list, night_number):
        #each attribute is declared here first, even if value is assigned later
        self.data = img_list #Read the list of images per nightnigh
        self.obs_night = night_number #each night is assigned a number, the first night is night 0
        self.image_data = None #the image data for the night
        self.headers = None #the night's headers
        self.wcs = None #world coordinate system transformation matrix
        self.readout_noise = None #detector readout noise
        self.aligned = None #aligned image data
        self.template = None #template frame
        self.is_reference = None #tf tag
        self.references = None #reference sources
        self.obs_filter = None #observation filter
        self.no_obs_images = None #number of images in night


    def initialize_frames(self):
        hdus = [image for image in self.data] #opens each fits image in a list using list comprehension
        self.image_data = [image['ALGN'].data for image in hdus] #pulls image data from file
        self.is_reference = self.image_data[0]
        self.headers = [image['ALGN'].header for image in hdus] #pulls header data from file.
        self.wcs = WCS(self.headers[0]) #gets wcs transformation, each night should only have one unique transformation
        self.readout_noise = self.headers[0]['RDNOISE'] #pulls readout noise from header. Readout noise is detector based and should be the same across nights taken with same equipment but do this just in case.
        self.mjd_times = []
        self.date_times = []
        self.date = []
        self.start_times = []
        for header in self.headers:
            self.mjd_times.append(header['MJD-OBS']) #Modified Julian Date (date+start time)
            self.date_times.append(header['DATE-OBS']) #Date/Start time (YYYY-MM-DD:HH:MM:SS in UTC)
            self.date.append(header['DAY-OBS']) #Date (YYYYMMDD)
            self.start_times.append(header['UTSTART']) #UTC start time HH:MM:SS

        self.template = np.median(self.image_data, axis = 0) #night template
        background = sep.Background(self.template) #sep background subtraction for source extraction
        self.references = sep.extract(self.template - background.back(),  background.globalrms*3, minarea =25, segmentation_map=False) #finds sources based on given parameters
        self.obs_filter = self.headers[0]['filter'][0] #observation filter.
        self.no_obs_images = len(self.image_data) #number of images in night.

    # def get_info(self): #function to grab night info, useful for debugging.
    #     print(f"path: {self.path}, night {self.obs_night}, n_frames: {len(self.image_data)}, n_aligned: {len(self.image_data)}, wcs: {self.wcs}, n_ref: {len(self.references)}, filter: {self.obs_filter}")

class Source: #initialize source object
    def __init__(self, source, count, WCS):
        self.position = pixel_to_skycoord(source['x'], source['y'], wcs= WCS).transform_to('icrs') #since pixel locations are inconsistent, store position as RA/DEC
        self.radius = (source['xmax'] - source['xmin'])/2 #source radius (size) provided by SEP
        self.source_id = count #identifying number
        self.is_reference = None #if star is reference
        self.ref_mag = None #SDSS magnitude, if available
        self.ref_mag_err = None #reference mag error
        self.inst_mags = [] #instrumnetal (our) magnitudes
        self.inst_mag_errs = [] #instrumental mag errors
        self.calibrated_mags = [] #calibrated magnitudes
        self.flagged = False #bad source flag. Will be flipped true if a source is not present in all observing nights or has negative aperture sum
        self.chi2 = None

    def query_source(self, night): #querry a source through the sdss database
        #we want the search to return ra, dec, mags, and mag error. region = False, returns first result of search.
        search = SDSS.query_crossid(self.position, fields = ['ra', 'dec', f'psfMag_{night.obs_filter}', f'psfMagErr_{Nights[0].obs_filter}'], radius = 15 * u.arcsec, region = False)
        if search:
            if search['type'] == "STAR": #want to make sure reference objects are Stars
                self.is_reference = True #lets us know star is a reference star
                self.ref_mag = search[f'psfMag_{night.obs_filter}'] #fill in mag and error fields:
                self.ref_mag_err = search[f'psfMagErr_{night.obs_filter}']

    def boundary_check(self, night): #checks if a star is within frame for a given night.
        source_xy = SkyCoord.to_pixel(self.position, wcs= night.wcs)
        if (night.headers[0]['NAXIS1'] - source_xy[0]) < 0 or source_xy[0] < 0 or (night.headers[0]['NAXIS2'] - source_xy[1]) < 0 or source_xy[1] < 0:
            self.flagged = True #if star is out of bounds, flags star as bad

    def aperture_photometry(self, img, nght): #does aperture photometry
        coords = SkyCoord.to_pixel(self.position, wcs = nght.wcs) #gets pixel values of source from RA DEC
        pcoords = PixCoord(coords[0], coords[1]) #another coord object needed for Regions
        radius_i = self.radius #inner aperture radius
        radius_o_0 = radius_i + 5 #inner annulus radius
        radius_o_1 = radius_o_0 + 5 #outer annulus radius
        source_circle = CirclePixelRegion(pcoords, radius_i).to_mask() #makes region of source shape
        source_aperture = source_circle.cutout(img) #gets data of source
        background_annulus = CircularAnnulus(coords, radius_o_0, radius_o_1) #makes annulus for background subtraction
        background_mean = ApertureStats(img, background_annulus).mean #takes mean of background annulus
        source_flux_pix = source_aperture-(source_circle*background_mean) #pixel wise background subtraction
        source_flux_total = np.sum(source_flux_pix) #total flux
        readout_sum_source = np.pi*(radius_i**2)*(nght.readout_noise**2)
        readout_sum_annulus = np.pi*((radius_o_1**2)-(radius_o_0**2))*(nght.readout_noise**2)
        delta_n = (readout_sum_source + source_flux_total + (((radius_i**2)/((radius_o_1**2)-(radius_o_0**2)))**2)*(readout_sum_annulus + aperture_photometry(img, background_annulus)['aperture_sum'][0]))**(1/2) #term needed for magnitude error calculation
        if source_flux_total < 0:
            print(self.source_id)
            self.flagged = True #flags source if aperture sum turns out to be negative
        else:
            instrumental_mag = -2.5*np.log10(source_flux_total) #magnitude
            instrumental_mag_error = 2.5*np.log10(np.e)*abs(delta_n/source_flux_total) #magntiude error
            self.inst_mags.append(instrumental_mag)
            self.inst_mag_errs.append(instrumental_mag_error)

    def add_calibrated_mag(self, mag):
        self.calibrated_mags.append(mag) #adds calibrated mag. For some reason, math comes out unexpectedly if calibration takes place in class function.

    def add_chi(self, chi):
        self.chi2 = chi

    def get_info(self): #prints out source info.
        print(f"path: {self.path}, night {self.obs_night}, n_frames: {len(self.image_data)}, n_aligned: {len(self.aligned_images)}, n_ref: {len(self.references)}, filter: {self.obs_filter}")


def lin_model(p, x): #define a standard linear model for ODR fitting. Part of calibration.
    return p[0] * x + p[1]

# def photometry_test(hduls):
#     start = time.perf_counter()
#     median_img = combine(hduls,'SUB')
#     median_img = median_img['PRIMARY'].data
#     median_img = median_img.astype(np.float64)
#     threshold = detect_threshold(median_img, nsigma=10.0,background=0.0)
#     segmentedimg = detect_sources(median_img, threshold=threshold,npixels=10)
#     sourcecatalog = SourceCatalog(median_img, segmentedimg)
#     print(' ')
#     end = time.perf_counter()
#     #for source in sourcecatalog:
#     #    print("Source Position: ", source.centroid)
#     print('Sources detected = {}. Time to detect sources = {} seconds'.format(len(sourcecatalog),end-start))
#     return (h for h in hduls)

def photometry(hduls, name, directory):
    Nights = [Night(h,i) for i,h in enumerate(hduls)]
    for night in Nights:
        night.initialize_frames() #see initialize_frames() in Night class definitions
        # night.get_info()

    Sources = [Source(source, count, Nights[0].wcs, Nights) for count, source in enumerate(Nights[0].references)] #initializes sources based off first night's list. This ensures proper source tracking
    for source in Sources:
        if source.flagged == False:
            source.query_source(Nights[0])  # see query_source() in source class definitions
    # print('len of sources is: {}'.format(len(Sources)))
    night_array = [] #this is to help organize plotting later.
    mag_thresh = 15 #magnitude threshold for calibrating sources.
    for night in Nights: #for each night, iterates through every source for each image.
        for image in night.image_data:
            for source in Sources:
                source.boundary_check(night) #see boundary_check() in source class definition
                if source.flagged == False:
                    source.aperture_photometry(image, night)  #see aperture_photometry() in source class definition
            night_array.append(night.obs_night)
   # print('Night array is {}'.format(night_array))
    slopes = []
    zeros = []
    counter = 0
    # slope_errs = []
    # zero_errs = []

    for night in Nights:
        for image in night.aligned_images:
            inst_mags = [source.inst_mags[counter] for source in Sources if source.is_reference == True and source.ref_mag < mag_thresh and source.flagged == False]
            inst_errs = [source.inst_mag_errs[counter] for source in Sources if source.is_reference == True and source.ref_mag < mag_thresh and source.flagged == False]
            ref_mags = [source.ref_mag[0] for source in Sources if source.is_reference == True and source.ref_mag < mag_thresh and source.flagged == False]
            ref_mag_errs = [source.ref_mag_err[0] for source in Sources if source.is_reference == True and source.ref_mag < mag_thresh and source.flagged == False]
            linear = odr.Model(lin_model)
            calibration_data = odr.Data(inst_mags, ref_mags, we=inst_errs, wd=ref_mag_errs)
            fit_params = odr.ODR(calibration_data, linear, beta0=[1.0, 23.5]).run()  # beta0 is initial guesses
            slopes.append(fit_params.beta[0])
            zeros.append(fit_params.beta[1])
            counter += 1

    #print(f"Slopes: {slopes[5]}")
    #print(f"Zero: {zeros[5]}")
    #print(f"Slope ERRs: {slope_errs[5]}")
    #print(f"Zero errs: {zero_errs[5]}")

    for source in Sources:
        if source.flagged != True:
            for i in range(0, len(slopes)):
                mag = (source.inst_mags[i]*slopes[i] + zeros[i])
                source.add_calibrated_mag(mag)
    
    color_arr = np.array(sorted(CSS4_COLORS, key=lambda c: tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(c)))))
    l = np.arange(0, len(slopes))

    median_curves = []
    for night in Nights:
        median_mags = []
        index_l = len(night.image_data)*night.obs_night
        index_h = len(night.image_data)*(1+night.obs_night)
        for source in Sources:
            if source.flagged != True:
                median_mags.append(source.calibrated_mags[index_l:index_h])
        median_curve = np.median(median_mags, axis = 0)/np.median(median_mags)
        median_curves.append(median_curve)
    med_curve = np.concatenate(median_curves)
    plt.scatter(l, med_curve, c = color_arr[93::5][np.array(night_array)])
    plt.gca().invert_yaxis()
    #PLT.SAVE

    for source in Sources:
        if source.flagged != True:
            if len(Nights) > 1:
                print('length is accurate')
                #print(source.get_info())
                night_array = np.array(night_array)
                source_mags = np.array(source.calibrated_mags/med_curve)
                source_errs = np.array(source.errors)
                median_mags = []
                median_errs = []
                mjd_times = []
                utc_times = []
                for night in range(np.max(night_array) + 1): # +1 avoids fencepost error in the night because range is exclusive.
                    print('Night number is: {}'.format(night))
                    #Magnitudes
                    index_array = np.where(night_array == night)[0]
                    night_mags = source_mags[index_array]
                    night_errs = source_errs[index_array]
                    median_mags.append(np.median(night_mags))
                    median_errs.append(np.median(night_errs))
                    #Times
                    mjd_array = np.array(Nights[night].mjd_times)
                    utc_array = np.array(Nights[night].date)
                    mjd_times.append(np.median(mjd_array))
                    print('time_array: {}'.format(mjd_times))
                    utc_times.append(utc_array[int(len(utc_array)/2)])
                # print(median_errs)
                avg_mag = np.average(source.calibrated_mags/med_curve, weights= source.weights)
                Chis = []
                for i, m in enumerate(source.calibrated_mags):
                    chi_i = (((m/med_curve[i]) - avg_mag)**2)/(source.errors[i]**2)
                    Chis.append(chi_i)
                dof = len(source.calibrated_mags) - 1
                chi_dof = np.sum(Chis)/dof
                source.add_chi(chi_dof)
                dof_string =  "%.2f" % chi_dof
                plt.figure(figsize=(10,7.5))
                plt.title(f"Source: {source.source_id}, Location: {source.position}, Chi2/Dof: {dof_string}")
                plt.errorbar(mjd_times, median_mags, yerr=median_errs, elinewidth=2, capsize=5, linestyle="", marker="o", color="black")
                plt.xticks(mjd_times, utc_times, rotation=45, ha='right')
                plt.xlabel("Observation Date (YYYYMMDD)")
                plt.ylabel(f"Magnitude, {Nights[0].obs_filter}-band")
                plt.gca().invert_yaxis()
                plt.tight_layout()
                plt.plot(mjd_times, np.ones(len(mjd_times)) * avg_mag, linestyle='--', color='black',
                         label="TRIPP Average Mag: {}".format("%.3f" % avg_mag))
                if source.is_reference:
                plt.plot(mjd_times, np.ones(len(mjd_times)) * source.ref_mag, linestyle='dashdot',
                         color=f"{Nights[0].obs_filter}", label="SDSS Mag: {}".format("%.3f" % source.ref_mag))                
                plt.legend()
                plt.savefig(directory + f"Source {source.source_id} Lightcurve.png", dpi = 1000)

            else: #single night
                # print(source.get_info())
                r = np.arange(0, len(source.calibrated_mags))
                avg_mag = np.average(source.calibrated_mags/med_curve, weights=source.weights)
                #avg_mag = np.mean(source.calibrated_mags/med_curve)
                Chis = []
                for i, m in enumerate(source.calibrated_mags):
                    chi_i = (((m / med_curve[i]) - avg_mag) ** 2) / (source.errors[i] ** 2)
                    Chis.append(chi_i)
                dof = len(source.calibrated_mags) - 1
                chi_dof = np.sum(Chis) / dof
                dof_string = "%.2f" % chi_dof
                plt.figure(figsize=(12, 9))
                plt.errorbar(r, source.calibrated_mags/med_curve, yerr=source.errors, elinewidth=1, capsize=2, markersize = 3, linestyle = 'none', marker = 'o', c = 'black')
                plt.plot(r, np.ones(len(r))*avg_mag, linestyle = '--', color = 'black', label = f"TRIPP Avg Mag:{avg_mag}")
                if source.is_reference:
                    plt.plot(r, np.ones(len(r))*source.ref_mag, linestyle = 'dashdot', color = f"{Nights[0].obs_filter}", label = "SDSS Mag: {}".format("%.3f" % source.ref_mag))
                x_locs= []
                times = []
                for i in range(len(Nights[0].image_data))[::10]:
                    x_locs.append(i)
                    times.append(Nights[0].start_times[i])
                plt.xlabel("Observation Start Time, UTC (HH:MM:SS)")
                plt.ylabel(f"Magnitude {Nights[0].obs_filter}-band")
                plt.title(f"Source: {source.source_id}, Location: {source.position}, Chi2/Dof: {chi_dof}")
                plt.gca().invert_yaxis()
                plt.xticks(x_locs, times, rotation = 45)
                plt.legend()
                os.chdir(directory)
                plt.savefig(directory + f"Source {source.source_id} Lightcurve.png", dpi = 1000)

    transient_candidates = []
    for source in Sources:
        if source.flagged != True:
            if source.chi2 > 25:
                transient_candidates.append(source.get_info())

    return (h for h in hduls)

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
    return photometry(hduls,name,directory)
