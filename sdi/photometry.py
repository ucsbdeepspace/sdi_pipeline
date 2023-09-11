import csv
import glob
import sep
import astropy.units as u
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import numpy as np
import astroalign as aa
from scipy import odr
import matplotlib.colors as mcolors
from astropy.io import fits
from astroquery.ipac.irsa import Irsa
from astroquery.sdss import SDSS
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy.coordinates import SkyCoord
from matplotlib.colors import LogNorm, TABLEAU_COLORS, CSS4_COLORS
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus, ApertureStats
from regions import CirclePixelRegion, PixCoord
from photutils import centroids
import multiprocessing as mp
from multiprocessing import shared_memory
from . import _cli as cli
import click
import os
from .combine import combine
from photutils.aperture import CircularAperture, aperture_photometry
from photutils import background as bck
from photutils.segmentation import detect_sources,detect_threshold,SourceCatalog
from photutils import detection
from astropy.stats import sigma_clipped_stats, SigmaClip
import time

class Night: #Each observing night is initialized as a Night object
    def __init__(self, img_list, night_number):
        #each attribute is declared here first, even if value is assigned later
        self.data = img_list #Read the list of images per night
        self.obs_night = night_number #each night is assigned a number, the first night is night 0
        self.image_data = None #the image data for the night
        self.headers = None #the night's headers
        self.wcs = None #world coordinate system transformation matrix
        self.readout_noise = None #detector readout noise
        self.aligned = None #aligned image data
        self.template = None #template frame
        self.references = None #reference sources
        self.obs_filter = None #observation filter
        self.no_obs_images = None #number of images in night


    def initialize_frames(self):
        hdus = [image for image in self.data] #opens each fits image in a list using list comprehension
        self.image_data = [image['ALGN'].data for image in hdus] #pulls image data from file
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
        self.aligned_images = self.image_data

        self.template = np.median(self.aligned_images, axis = 0) #night template
        background = sep.Background(self.template) #sep background subtraction for source extraction
        self.references = sep.extract(self.template - background.back(),  background.globalrms*3, minarea =25, segmentation_map=False) #finds sources based on given parameters

        self.obs_filter = self.headers[0]['filter'][0] #observation filter.
        self.no_obs_images = len(self.aligned_images) #number of images in night.

    def get_info(self): #function to grab night info, useful for debugging.
        print(f"path: {self.path}, night {self.obs_night}, n_frames: {len(self.image_data)}, n_aligned: {len(self.aligned_images)}, wcs: {self.wcs}, n_ref: {len(self.references)}, filter: {self.obs_filter}")

class Source: #initialize source object
    def __init__(self, source, count, WCS, Nights):
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
        self.weights = []
        self.errors = []
        self.chi2 = None
        self.Nights = Nights


    def query_source(self): #querry a source through the sdss database
        #we want the search to return ra, dec, mags, and mag error. region = False, returns first result of search.
        Nights = self.Nights
        search = SDSS.query_crossid(self.position, fields = ['ra', 'dec', f'psfMag_{Nights[0].obs_filter}', f'psfMagErr_{Nights[0].obs_filter}'], radius = 15 * u.arcsec, region = False)
        if search:
            if search['type'] == "STAR": #want to make sure reference objects are Stars
                self.is_reference = True #lets us know star is a reference star
                self.ref_mag = search[f'psfMag_{Nights[0].obs_filter}'] #fill in mag and error fields:
                self.ref_mag_err = search[f'psfMagErr_{Nights[0].obs_filter}']

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

        source_circle = CirclePixelRegion(pcoords, radius_i) #makes region of source shape
        source_circle_mask = source_circle.to_mask()
        source_aperture = source_circle_mask.cutout(img) #gets data of source
        source_sum_unsub = np.sum(source_aperture)

        background_annulus = CircularAnnulus(coords, radius_o_0, radius_o_1) #makes annulus for background subtraction
        #background_mean = ApertureStats(img, background_annulus).mean #takes mean of background annulus
        background_sum = aperture_photometry(img, background_annulus)['aperture_sum'][0]

        #source_flux_pix = source_aperture-((source_circle.area/background_annulus.area)*background_sum*source_circle_mask) #pixel wise background subtraction
        source_flux_total = np.sum(source_aperture) - (source_circle.area/background_annulus.area)*background_sum  #total flux



        readout_sum_source = source_circle.area*(nght.readout_noise**2)
        readout_sum_annulus = background_annulus.area*(nght.readout_noise**2)

        delta_n = (readout_sum_source + source_flux_total + ((source_circle.area/background_annulus.area)**2)*(readout_sum_annulus + background_sum))**(1/2)

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

    def add_error(self, err):
        self.errors.append(err)
        self.weights.append(1/(err**2))

    def get_info(self): #prints out source info.
        Nights = self.Nights
        print(f"ra_dec: {self.position}, Night_0_xy: {SkyCoord.to_pixel(self.position, wcs = Nights[0].wcs)} rad: {self.radius}, ref_status: {self.is_reference}, ref_mag: {self.ref_mag}, inst_mag_avg:{np.mean(self.inst_mags)}, cal_mag_avg: {np.mean(self.calibrated_mags)}, flagged: {self.flagged}, ID: {self.source_id}, Chi2: {self.chi2}")

    def __iter__(self): #for writing out csv files
        return iter([self.position, self.is_reference, self.ref_mag, self.chi2, self.flagged, self.source_id, self.calibrated_mags])


def lin_model(p, x): #define a standard linear model for ODR fitting. Part of calibration.
    return p[0] * x + p[1]

def ODR(x_data, y_data):
    x_bar = np.mean(x_data)
    y_bar = np.mean(y_data)

    s_xx = 1/len(x_data) * np.sum((x_data - x_bar)**2)
    s_yy = 1/len(y_data) * np.sum((y_data - y_bar)**2)
    s_xy = 1/len(x_data) * np.sum((x_data - x_bar) * (y_data - y_bar))

    b_0 = (s_yy - s_xx + np.sqrt((s_yy - s_xx)**2 + 4*s_xy**2))/(2 * s_xy)
    b_1 = y_bar - b_0 * x_bar

    return [b_0, b_1]

def photometry_test(hduls):
    start = time.perf_counter()
    median_img = combine(hduls,'SUB')
    median_img = median_img['PRIMARY'].data
    median_img = median_img.astype(np.float64)
    threshold = detect_threshold(median_img, nsigma=10.0,background=0.0)
    segmentedimg = detect_sources(median_img, threshold=threshold,npixels=10)
    sourcecatalog = SourceCatalog(median_img, segmentedimg)
    print(' ')
    end = time.perf_counter()
    #for source in sourcecatalog:
    #    print("Source Position: ", source.centroid)
    print('Sources detected = {}. Time to detect sources = {} seconds'.format(len(sourcecatalog),end-start))
    return (h for h in hduls)

def photometry(hduls, name, directory):

    Nights = [Night(h,i) for i,h in enumerate(hduls)]
    for night in Nights:
        night.initialize_frames() #see initialize_frames() in Night class definitions

    Sources = [Source(source, count, Nights[0].wcs, Nights) for count, source in enumerate(Nights[0].references)] #initializes sources based off first night's list. This ensures proper source tracking
    for source in Sources:
        source.query_source() #see query_source() in source class definitions
    print('len of sources is: {}'.format(len(Sources)))
    night_array = [] #this is to help organize plotting later.
    mag_thresh = 15 #magnitude threshold for calibrating sources.
    for night in Nights: #for each night, iterates through every source for each image.
        for image in night.aligned_images:
            for source in Sources:
                source.boundary_check(night) #see boundary_check() in source class definition
                if source.flagged == False:
                    source.aperture_photometry(image, night)  #see aperture_photometry() in source class definition
            night_array.append(night.obs_night)
   # print('Night array is {}'.format(night_array))
    slopes = []
    zeros = []
    slope_errs = []
    zero_errs = []

    counter = 0

    for night in Nights:
        for image in night.aligned_images:
            #print('image: '.format(image))
            instrumental_magnitudes = [s.inst_mags[counter] for s in Sources if s.is_reference == True and s.ref_mag < mag_thresh and s.flagged != True]
            print('inst magnitude is {}'.format(instrumental_magnitudes))
            reference_magnitudes = [s.ref_mag[0] for s in Sources if s.is_reference == True and s.ref_mag < mag_thresh and s.flagged != True]
            print('reference magnitude is: {}'.format(reference_magnitudes))
            jk_params = np.zeros((len(instrumental_magnitudes), 2))
            for i in range(len(instrumental_magnitudes)):
                x_sample = np.append(instrumental_magnitudes[:i], instrumental_magnitudes[i+1:])
                y_sample = np.append(reference_magnitudes[:i], reference_magnitudes[i+1:])
                jk_params[i] = ODR(x_sample, y_sample)
                #print('jk params:{}'.format(jk_params[i]))

            mean_params = np.mean(jk_params, axis = 0)
            sig_params = np.std(jk_params, axis = 0)

            slopes.append(mean_params[0])
            zeros.append(mean_params[1])
            slope_errs.append(sig_params[0])
            zero_errs.append(sig_params[1])
            print('slopes are now: {}'.format(slopes))
            counter += 1

    #print(f"Slopes: {slopes[5]}")
    #print(f"Zero: {zeros[5]}")
    #print(f"Slope ERRs: {slope_errs[5]}")
    #print(f"Zero errs: {zero_errs[5]}")

    for source in Sources:
        if source.flagged != True:
            for i in range(0, len(slopes)):
                mag = (source.inst_mags[i]*slopes[i] + zeros[i])
                # final_err = np.sqrt((slopes[i]* source.inst_mags[i])**2 * ((slope_errs[i]/slopes[i])**2 + (source.inst_mag_errs[i]/source.inst_mags[i])**2) + zero_errs[i]**2)
                final_err = source.inst_mag_errs[i]
                source.add_calibrated_mag(mag)
                source.add_error(final_err)
    
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
                print(median_errs)
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
                plt.plot(mjd_times, np.ones(len(mjd_times))*avg_mag, linestyle = '--', color = 'black', label = "TRIPP Average Mag: {}".format("%.3f" % avg_mag))
                if source.is_reference:
                    plt.plot(mjd_times, np.ones(len(mjd_times))*source.ref_mag, linestyle = 'dashdot', color = f"{Nights[0].obs_filter}", label = "SDSS Mag: {}".format("%.3f" % source.ref_mag))
                plt.legend()
                os.chdir(directory)
                plt.savefig(f"Source {source.source_id} Lightcurve.png", dpi = 1000)

            else: #single night
                print(source.get_info())
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
                plt.savefig(f"Source {source.source_id} Lightcurve.png", dpi = 1000)
    os.chdir(directory)
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
