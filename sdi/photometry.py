import glob
import click
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
from matplotlib.colors import LogNorm, TABLEAU_COLORS
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus, ApertureStats
from regions import CirclePixelRegion, PixCoord
from photutils import centroids
from . import _cli as cli

def photometry(files_dir, out_dir , filter , plot_mode, median_divide): #barebones, can add more stuff later
    files = sorted(glob.glob(files_dir + '/*'))  # filepath, this is a directory containing directories.
    outstr= out_dir
    print(outstr)
    data = []
    for i, dir in enumerate(files):
        set = dict(images=[], header=None, template=None, refs=None, WCS=None,run=i)  # dictionary object used to sort input files
        images = sorted(glob.glob(dir + '/*'))  # formatting
        print(dir)
        hdus = [fits.open(i) for i in images]  # opens fits files
        frames = [h[1].data for h in hdus]  # image data
        header = (hdus[0])[1].header  # the header for each run is the header for the first image, you don't need each individual header.
        try:  # some datsets won't align, this is not ideal but we can skip alignment.
            aligned = [aa.register(i, frames[0])[0] for i in frames[0:]]  # use astroalign for the time being, want to use multiprcoessed eventually
        except:
            aligned = frames
            print("DID NOT ALIGN")

        template = np.median(aligned, axis=0)  # creates median value template
        w = WCS(hdus[0][1].header)  # WCS matrix object, used to transform pixel values to RA-DEC coordinates

        # Fills out set object
        set['images'] = aligned
        set['header'] = header
        set['template'] = template
        set['WCS'] = w
        bkg_phot = sep.Background(template)  # background subtract for source extraction.
        set['refs'] = sep.extract(template - bkg_phot.back(), bkg_phot.globalrms * 3, minarea=25, segmentation_map=False)  # find sources in image
        # this ensures that the source list is the same for all sets
        set['run'] = i
        print(len(set['refs']))  # check
        data.append(set)

    sources = []
    n_ref = 0  # check counter
      # observation filter
    for c, src in enumerate(data[0]['refs']):  # looking at the sources picked out in the first run so that everything is with respect to one source list.
        x = src['x']
        y = src['y']
        coord = pixel_to_skycoord(x, y, data[0]['WCS']).transform_to('icrs')  # gives wcs transformation for pixel coordinates
        search = SDSS.query_crossid(coord,fields=['ra', 'dec', 'psfMag_{}'.format(filter), 'psfMagErr_{}'.format(filter)],radius=15 * u.arcsec, region=False)  # narrow field cone search to find source based on ra, dec.
        radius = (src['xmax'] - src['xmin']) / 2
        ref = dict(ra_dec=coord, source_id=c, rad=radius, ref_mag=None, ref_mag_err=None,border=False)  # because the pixel value isn't static across observing runs, don't save the x-y values just ra dec and convert later.
        for i in range(0, len(data)):  # creates mag slots for each observing run
            ref['calibrated_mags_{}'.format(i)] = []  # even though the calibrated mags should line up (and not need to be separated) this is to help visualize day changes later on
            ref['instrumental_mags_{}'.format(i)] = []
            ref['inst_mag_errs_{}'.format(i)] = []
        if search:  # if SDSS query returned results, continue
            if search['type'] == 'STAR':  # filters search results by sources that of type star.
                ref['ref_mag'] = search['psfMag_{}'.format(filter)]  # add reference mag
                ref['ref_mag_err'] = search['psfMagErr_{}'.format(filter)]  # add SDSS error
                n_ref += 1
        sources.append(ref)

    for set in data:
        print(set['run'])
        for i, image in enumerate(set['images']):
            print(i)
            N_r = set['header']['RDNOISE']  # readout noise
            for source in sources:
                coords = SkyCoord.to_pixel(source['ra_dec'], wcs=set['WCS'])  # gets pixel values of source from RA DEC
                # print(coords[0], coords[1])
                if (set['header']['NAXIS1'] - coords[0]) < 0 or coords[0] < 0 or (set['header']['NAXIS2'] - coords[1]) < 0 or coords[1] < 0:  # checks to see if a source exceeds image boundaries for later image sets.
                    source['border'] = True
                    print(source['source_id'], coords[0], coords[1])
                    continue

                pcoords = PixCoord(coords[0], coords[1])  # another coord object needed for Regions
                radius_i = source['rad']
                radius_o_0 = radius_i + 5  # inner annulus radius
                radius_o_1 = radius_o_0 + 5  # outer annulus radius

                source_circle = CirclePixelRegion(pcoords, radius_i).to_mask()  # makes region of source shape
                source_aperture = source_circle.cutout(image)  # gets data of source

                background_annulus = CircularAnnulus(coords, radius_o_0, radius_o_1)
                background_mean = ApertureStats(image, background_annulus).mean

                source_flux_pix = source_aperture - (source_circle * background_mean)  # pixel wise background subtraction
                source_flux_total = np.sum(source_flux_pix)  # total flux

                readout_sum_square = np.sum(source_circle * np.float64(N_r ** 2))  # applies square readout noise to source array shape, then adds. Gives sum of square readout noise over back subtracted source.
                delta_n = (readout_sum_square + source_flux_total + (((radius_i ** 2) / ((radius_o_1 ** 2) - (radius_o_0 ** 2))) ** 2) * (readout_sum_square + aperture_photometry(image, background_annulus)['aperture_sum'][0])) ** (1 / 2)  # this is the stuff for SNR
                if source_flux_total <= 0:
                    inst_mag = -2.5 * np.log10(
                        abs(source_flux_total))  # For now, the case where the background is oversubtracted from LCO is handled in this way but this is probably not the correct way to do this.
                    delta_m = 2.5 * np.log10(np.e) * abs(delta_n / source_flux_total)  # magnitude error

                else:
                    inst_mag = -2.5 * np.log10(source_flux_total)
                    delta_m = 2.5 * np.log10(np.e) * abs(delta_n / source_flux_total)
                source['instrumental_mags_{}'.format(set['run'])].append(inst_mag)
                source['inst_mag_errs_{}'.format(set['run'])].append(delta_m)

    # NEED TO PROPAGATE ERROR FROM THIS SECTION

    for set in data:
        res = []
        mag_thresh = 15  # magnitude threshold for picking reference stars

        # criteria for reference stars are that they have an SDSS reference mag, are brighter than the mag threshold and are not outside of the image
        inst_mags = [np.mean(source['instrumental_mags_{}'.format(set['run'])]) for source in sources if source['ref_mag'] != None and source['ref_mag'] < mag_thresh and source['border'] == False]
        sky_mags = [source['ref_mag'][0] for source in sources if source['ref_mag'] != None and source['ref_mag'] < mag_thresh and source['border'] == False]
        # Makes linear model for calibration:
        # This is the first round of modeling, with outliers.
        print(len(inst_mags), len(sky_mags))
        p0 = np.polyfit(inst_mags, sky_mags, deg=1)  # linear fit for instrumental to SDSS mags
        x = np.arange(-15, 0)  # plotting fit line
        y = p0[0] * x + p0[1]  # fit line
        plt.plot(x, y, color='b', label="Model With Outliers, m = {}, b = {}".format("%.4f" % p0[0], "%.4f" % p0[1]))
        diffs = [s['ref_mag'][0] - (np.mean(s['instrumental_mags_{}'.format(set['run'])]) * p0[0] + p0[1]) for s in sources if s['ref_mag'] != None and s['ref_mag'] < mag_thresh and s['border'] == False]
        stdv = np.std(diffs)
        inst_mags_final = []
        sky_mags_final = []
        outlier_inst = []
        outlier_sky = []

        for diff in diffs:  # rudementary sigma clipping to remove outliers from calibration model.
            if diff < stdv:
                i = diffs.index(diff)
                inst_mags_final.append(inst_mags[i])
                sky_mags_final.append(sky_mags[i])
            else:
                i = diffs.index(diff)
                outlier_inst.append(inst_mags[i])
                outlier_sky.append(sky_mags[i])
        p1 = np.polyfit(inst_mags_final, sky_mags_final, deg=1)  # recalculates calibration model without outliers.
        # p2 = np.polyfit(inst_mags_final, sky_mags_final, deg = 0)
        # print(p2[0])
        print("first try: {}".format(p0))  # prints slopes of each model. In theory, they should come out to around 1.
        print("second try: {}".format(p1))

        plt.scatter(outlier_inst, outlier_sky, color='b', s=4, label="Outliers")
        plt.scatter(inst_mags_final, sky_mags_final, color='r', s=4, label="Kept")
        plt.plot(x, [i * p1[0] + p1[1] for i in x], color='r',
                 label="Model Without Outliers, m = {}, b = {}".format("%.4f" % p1[0], "%.4f" % p1[1]))
        # plt.plot(x, [i+ p1[1] for i in x], color = 'g', label = "unity")
        plt.xlabel("Instrumental Magnitude SDSS g-band")
        plt.ylabel("SDSS Reference Magnitude g-band")
        plt.title("Instrumental vs Reference Magnitude")
        plt.legend()
        #add dynamic plot saving back in here.
        plt.show()

        # add calibrated mags to sources:
        for source in sources:
            vals = []
            for val in source['instrumental_mags_{}'.format(set['run'])]:
                if abs(1 - p1[0]) < abs(1 - p0[0]):
                    cal = np.float64(val * p1[0] + p1[1])
                    vals.append(cal)
                else:
                    cal = np.float64(val * p0[0] + p0[1])
                    vals.append(cal)
            source['calibrated_mags_{}'.format(set['run'])] = vals  # probably a cleaner way to do this part but was having issue where calibrated magnitudes were being added to dict as individual arrays



    for set in data: #this part is entirely optional, need to nest in if statement.
        print(set['run'])
        differences = [s['ref_mag'] - np.mean(s['calibrated_mags_{}'.format(set['run'])]) for s in sources if s['ref_mag'] != None and s['border'] == False]
        mean_diff = np.mean(differences)
        std_diff = np.std(differences)
        med_diff = np.median(differences)
        print(mean_diff, med_diff, std_diff)

        X = [np.mean(s['calibrated_mags_{}'.format(set['run'])]) for s in sources if s['ref_mag'] != None and s['ref_mag'] < mag_thresh and s['border'] == False]
        Y = [s['ref_mag'] for s in sources if s['ref_mag'] != None and s['ref_mag'] < mag_thresh and s['border'] == False]
        p3 = np.polyfit(X, Y, deg=1)
        x_vals = np.arange(9.5, 18.5)
        plt.plot(x_vals, [x * p3[0] + p3[1] for x in x_vals], color='b', label='m = {}, b = {}'.format("%.4f" % p3[0], "%.4f" % p3[1]))

        print(p3)

        for s in sources:
            if s['ref_mag'] != None and abs(s['ref_mag'] - np.mean(s['calibrated_mags_{}'.format(set['run'])])) > 3 * std_diff:
                plt.errorbar(np.mean(s['calibrated_mags_{}'.format(set['run'])]), s['ref_mag'], xerr=np.median(s['inst_mag_errs_{}'.format(set['run'])]), yerr=s['ref_mag_err'],linestyle='none', marker='o', markersize=2.5, color='b')
                print(s['source_id'])
            else:
                plt.errorbar(np.mean(s['calibrated_mags_{}'.format(set['run'])]), s['ref_mag'], xerr=np.median(s['inst_mag_errs_{}'.format(set['run'])]), yerr=s['ref_mag_err'], linestyle='none', marker='o', markersize=2.5, color='r')

        plt.title("Calibrated vs Reference Magnitude")
        plt.xlabel("TRIPP Calibrated Magnitude SDSS g-band")
        plt.ylabel("SDSS Reference Magnitude g-band")
        plt.legend()
        #add dynamic plot saving back in.
        plt.show()


    ######
    runs = np.arange(0, len(data))
    colors = []
    for run in runs:  # sets up color scheme for plotting.
        i = 0
        for val in sources[0]['calibrated_mags_{}'.format(run)]:
            i += 1
            colors.append(run)
        print(i)
    c = np.array(colors)
    colors_options = np.array(['#3B5BA5', '#E87A5D', '#F3B941', '#f00a42', '#6F9BA4', '#9D9EA2', "#C5E86C", "#B4B5DF"])  # HEXIDECIMAL COLORS need to add a ton of colors to this even if they don't all get used. Could set it up do randomly generate hex colors but that will be inconsistent and kinda look like shit.
    #######
    mode = plot_mode
    med_dev = median_divide
    if med_dev == True:  # only run this for single night observing runs. Make part of each night?
        median_curves = []
        for i, run in enumerate(runs):
            median_mags = []
            for source in sources:
                if source['border'] == False:
                    median_mags.append([source[f'calibrated_mags_{i}']])
            median_curve = np.median(median_mags, axis=0) / np.median(median_mags)  # normalized median light curve
            median_curves.append(median_curve)
        curve = np.concatenate(median_curves)
    med_curve_final = np.concatenate(curve)
    p = np.arange(0, len(med_curve_final))
    plt.scatter(p, med_curve_final)
    plt.show()
    #######
    transient_candidates = []
    # plots light curves.
    for source in sources:
        x, y = SkyCoord.to_pixel(source['ra_dec'], wcs=data[0]['WCS'])
        tick_locs = []
        if mode == 0:
            MAGS = np.concatenate([source['calibrated_mags_{}'.format(i)] for i in runs if source['border'] == False])
            ERRS = np.concatenate([source['inst_mag_errs_{}'.format(i)] for i in runs if source['border'] == False])
            r = np.arange(0, len(MAGS), 1)
            set_counter = 0
            for i, val in enumerate(MAGS):
                if len(data) > 1 and i % len(data[set_counter]) == 0:  # THIS WILL ONLY WORK IF ALL SETS HAVE 25 FRAMES. Need to make it dynamic in case different nights don't have the same number of images.
                    plt.errorbar(i, val, yerr=ERRS[i], linestyle='none', marker='o', markersize=4, c=colors_options[c][i])
                    tick_locs.append(i)
                    set_counter += 1
                else:
                    plt.errorbar(i, val, yerr=ERRS[i], linestyle='none', marker='o', markersize=4,c=colors_options[c][i])
            if len(data) > 1:
                plt.xticks(tick_locs, [data[i]['header']['DAY-OBS'] for i in runs], rotation=20)

        if mode == 1:
            MAGS = [np.median(source['calibrated_mags_{}'.format(i)]) for i in runs if source['border'] == False]
            ERRS = [np.median(source['inst_mag_errs_{}'.format(i)]) for i in runs if source['border'] == False]
            r = np.arange(0, len(MAGS), 1)
            for i, mag in enumerate(MAGS):
                plt.errorbar(i, mag, yerr=ERRS[i], linestyle='none', marker='o', markersize=4, c=colors_options[r][i])
                tick_locs.append(i)
            plt.xticks(tick_locs, [data[i]['header']['DAY-OBS'] for i in runs], rotation=20)

        Chis = []
        avg_mag = np.mean(MAGS)
        for i, m in enumerate(MAGS):
            chi_i = ((m - avg_mag) ** 2) / (ERRS[i] ** 2)
            Chis.append(chi_i)
        dof = len(MAGS) - 1
        chi_dof = np.sum(Chis) / dof

        # transient/variable detection:
        dev = np.std(MAGS)
        print(dev)

        plt.title("X, Y: {}, {}; RA_DEC: {}; CHI2: {}, ID: {}".format("%.2f" % x, "%.2f" % y, source['ra_dec'],"%.2f" % chi_dof, source['source_id']))
        plt.plot(r, np.ones(len(r)) * avg_mag, label="TRIPP AVG MAG", linestyle='--', color='b')
        if source['ref_mag'] != None:
            if source['ref_mag'] < 16:
                plt.plot(r, np.ones(len(r)) * source['ref_mag'], linestyle='--', color='r', label="SDSS MAG [REF]")
            if source['ref_mag'] >= 16:
                plt.plot(r, np.ones(len(r)) * source['ref_mag'], linestyle='--', color='g', label="SDSS MAG [NOT REF]")
        plt.ylabel("TRIPP Magnitude, SDSS-{}'".format(filter))
        plt.xlabel("Observation Date YYYY/MM/DD")
        plt.legend()
        #dynamic plot saving
        plt.gca().invert_yaxis()
        #print("Plotted curve")
        plt.show()

        plt.title("Source Number: {}, Position: {}, {}".format(source['source_id'], "%.2f" % x, "%.2f" % y))
        circle0 = plt.Circle((x, y), source['rad'], color='r', fill=False)
        circle1 = plt.Circle((x, y), source['rad'] + 5, color='b', fill=False)
        circle2 = plt.Circle((x, y), source['rad'] + 5, color='g', fill=False)
        ax = plt.gca()
        ax.add_patch(circle0)
        ax.add_patch(circle1)
        ax.add_patch(circle2)
        plt.imshow(data[0]['template'], cmap='gray', norm=LogNorm(vmin=1, vmax=200), origin='lower')
        #dynamic plot saving here
        #print("plotted location")
        plt.show()
    ######

@cli.cli.command("photometry")
@click.option('-d', '--files_dir', type=str, help="Specify path to directory of fits files.")
@click.option('-D', "--out_dir", type=str, help = "Specify path to directory for output files", required = False)
@click.option('-f', "--filter", type = str, default = 'r', help = "Specify imaging filter for data set using SDSS color scheme", required = True)
@click.option('-p', "--plot_mode", type=int, default = 0 ,help ="Specify plotting mode", required = False)
@click.option('-m', "--median_divide", type=bool, default = True, help ="Divides by a median light curve to reduce systemic trends from light curves", required = False)
#@cli.operator

def photometry_cmd(files_dir, out_dir = "None", filter = 'r', plot_mode =0, median_divide = True):
    """
    does the photometry
    """
    return photometry(files_dir, out_dir, filter, plot_mode, median_divide)