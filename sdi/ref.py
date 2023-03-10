"""
ref -- this module matches sources to reference stars
HISTORY
    Created 2021-04-24
        Varun Iyer <varun_iyer@ucsb.edu>
"""
# general imports

import warnings
import click
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from . import _cli as cli

# define specific columns so we don't get dtype issues from the chaff
COLUMNS = ["source_id", "ra", "ra_error", "dec", "dec_error",
           "phot_g_mean_flux", "phot_g_mean_flux_error", "phot_g_mean_mag",
           "phot_rp_mean_flux", "phot_rp_mean_flux_error", "phot_rp_mean_mag",
           "phot_bp_mean_flux", "phot_bp_mean_flux_error", "phot_bp_mean_mag"]

def _in_cone(coord: SkyCoord, cone_center: SkyCoord, cone_radius: u.deg):
    """
    Checks if SkyCoord coord is in the cone described by conecenter and
    cone_radius
    """
    d = (coord.ra - cone_center.ra) ** 2 + (coord.dec - cone_center.dec) ** 2
    return d < (cone_radius ** 2)

def ref(hduls, read_ext="CAT", write_ext="REF"):
    """
    add information about remote reference stars to a 'REF' BinTableHDU
    \b
    Reference adds a new BinTableHDU entitled 'REF' which contains retrieved
    information about reference stars and their associations with the sources in
    'CAT'.  Each index of 'CAT' will correspond to the same index in 'REF'; e.g.
    reference star info associated with hduls[0]['CAT'].data[0] will be in
    hduls[0]['REF'].data[0]. If there is no associated reference information,
    any(hduls[0]['REF'].data[0]) will be False.
    :param hdul: A collection or generator of HDUL
    :param read_ext: the HDU extname to read source information from.
        Must include 'ra' and 'dec' fields.
    :param write_ext: the HDU extname to write reference information from.
    """
    # This import is here because the GAIA import slows down `import sdi`
    # substantially; we don't want to import it unless we need it
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = -1

    cached_table = np.array([])

    # Conduct a cone search that covers all possible sources in the image with one query
    # Use the header info in the first image for cone search size
    hduls = [h for h in hduls]
    hdul_init = hduls[0]
    cone_radius = max([hdul_init['ALGN'].header['NAXIS1'], hdul_init['ALGN'].header['NAXIS2']])*hdul_init['ALGN'].header['PIXSCALE']/3600
    cone_radius = u.Quantity(cone_radius, u.deg)
    # Make the query
    centercoord = SkyCoord(hdul_init['ALGN'].header['RA'], hdul_init['ALGN'].header['DEC'], unit="deg")
    data = Gaia.cone_search_async(centercoord, cone_radius, columns=COLUMNS, output_format="csv").get_results()
    data = data.as_array()
    cached_table = data.data
    print(len(cached_table))
    
    # we need this to track blanks till we know the dtype
    initial_empty = 0
    
    for hdul in hduls:
    # An adaptive method of obtaining the threshold value
        threshold = max(hdul[read_ext].data["a"])*hdul['ALGN'].header['PIXSCALE']/3600
        threshold = u.Quantity(threshold, u.deg)
        ra = hdul[read_ext].data["RA"]
        dec = hdul[read_ext].data["DEC"]
        output_table = np.array([])
        coordinates = SkyCoord(ra, dec, unit="deg")

        for coord in coordinates:
            ########### Query an area if we have not done so already ###########
            # Check to see if we've queried the area
            # Disable this whole section for now
            '''
            if not any((_in_cone(coord, query, cone_radius - 2 * threshold) \
                        for query in queried_coords)):
                # we have never queried the area. Do a GAIA cone search
                data = Gaia.cone_search_async(coord, cone_radius, columns=COLUMNS,
                                              output_format="csv").get_results()
                data = data.as_array()
                # add the cache table to the data
                if len(cached_table):
                    cached_table = np.hstack((cached_table, data.data))
                else:
                    cached_table = data.data
                for d in data:
                    # construct Coord objects for the new data
                    cached_coords.append(SkyCoord(d["ra"], d["dec"], unit=(u.deg, u.deg)))
                # note that we have now queried this area
                queried_coords.append(coord)
            '''

            ########### Look through our cache for matches #####################
            appended = False
            for ct, cs in zip(cached_table, coordinates):
                # look through the cache to find a match
                if _in_cone(cs, coord, threshold):
                    # if we find a match, copy it to the output table
                    if len(output_table):
                        output_table = np.hstack((output_table, np.copy(ct)))
                    else:
                        output_table = np.copy(ct)
                        output_table = np.hstack((np.empty(shape=initial_empty,
                                                           dtype=output_table.dtype), output_table))
                    appended = True
                    break
                else:
                    pass

            ########### Add a blank if we didn't find anything #################
            if not appended:
                if len(output_table):
                    # If we do not find one cached, then add a blank
                    blank = np.empty(shape=0, dtype=output_table.dtype)
                    output_table = np.hstack((output_table, blank))
                else:
                    initial_empty += 1

        ########## After going through all sources, add an HDU #################
        extname = write_ext
        header = fits.Header([fits.Card("HISTORY", "From the GAIA remote db")])
        # replace nan values with 0.0
        for i, elm in enumerate(output_table):
            for j, val in enumerate(elm):
                if np.isnan(val):
                    elm[j] = 0.0
        # only append the hdul if output_table is not empty
        if len(output_table):
            hdul.append(fits.BinTableHDU(data=output_table, header=header, name=extname))
        else:
            warnings.warn(f"empty reference table created, no stars found in the database corresponding to {hdul}")
        yield hdul
    return

@cli.cli.command("ref")
@click.option("-r", "--read-ext", default="CAT", help="The HDU to match")
@click.option("-w", "--write-ext", default="REF", help="The HDU to load ref into")
@cli.operator
def ref_cmd(hduls, read_ext="CAT", write_ext="REF"):
    """
    add information about remote reference stars to a 'REF' BinTableHDU
    \b
    Reference adds a new BinTableHDU entitled 'REF' which contains retrieved
    information about reference stars and their associations with the sources in
    'CAT'.  Each index of 'CAT' will correspond to the same index in 'REF'; e.g.
    reference star info associated with hduls[0]['CAT'].data[0] will be in
    hduls[0]['REF'].data[0]. If there is no associated reference information,
    any(hduls[0]['REF'].data[0]) will be False.
    :param hdul: A collection or generator of HDUL
    :param read_ext: the HDU extname to read source information from. Must include 'ra' and 'dec' fields.
    :param write_ext: the HDU extname to write reference information from.
    """
    try:
        read_ext = int(read_ext)
    except ValueError:
        pass
    try:
        write_ext = int(write_ext)
    except ValueError:
        pass
    return ref(hduls, read_ext, write_ext)
