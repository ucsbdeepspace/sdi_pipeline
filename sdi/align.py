"""
align -- this module aligns sets of astronomical data
HISTORY
    Created/Substantially refactored on 2019-09-01
        Varun Iyer <varun_iyer@ucsb.edu>
    Refactored SNR section on 2022-2-22
        Jasper Webb <jlw@ucsb.edu>
"""
import click
import numpy as np
import multiprocessing as mp
import multiprocessing.shared_memory as shared_memory
import os
import copy
import astroalign
#from .snr_function import snr
from photutils.aperture import CircularAperture, aperture_photometry
from photutils import background as bck
from photutils.segmentation import detect_sources,detect_threshold,SourceCatalog
from photutils import detection
from astropy.stats import sigma_clipped_stats, SigmaClip
import datetime
# from _scripts import snr
# from astropy.io.fits import CompImageHDU
# from sources import Source
from . import _cli as cli
import faulthandler

def multiprocessed_snr(read_ext, hduls, index, shm_name):
    hdul = hduls[index]
    try:
        data = hdul[read_ext].data
    except KeyError:
        data = hdul["PRIMARY"].data
    # identify background rms
    boxsize = (data.shape)
    bkg_estimator = bck.MedianBackground()
    sigma_clip = SigmaClip(sigma=4.)
    bkg = bck.Background2D(data, boxsize, bkg_estimator = bkg_estimator)
    bkg_mean_rms = np.mean(bkg.background_rms)

    data -= bkg.background    # set threshold and detect sources, threshold 5*std above background
    threshold = detect_threshold(data, nsigma=5.0,background=0.0)
    segmentedimg = detect_sources(data, threshold=threshold,npixels=10)
    sourcecatalog = SourceCatalog(data, segmentedimg)

    source_max_values = sourcecatalog.max_value
    avg_source_max_values = np.mean(source_max_values)
    existing_shm = shared_memory.SharedMemory(name = shm_name)
    snr_arr=np.ndarray(len(hduls), dtype = np.float64, buffer = existing_shm.buf)
    # calculate signal to noise ratio
    signal = avg_source_max_values
    noise = bkg_mean_rms
    sig_to_noise = (signal)/(noise)
    snr_arr[index] = sig_to_noise
    existing_shm.close()

def snr(shm_name, hduls, read_ext="SCI"):
    """
    calculates the signal-to-noise ratio of a fits file
    """
    processes = []
    for i in range(len(hduls)):
        p = mp.Process(target=multiprocessed_snr, args = (read_ext, hduls, i, shm_name))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()

def multiprocessed_align(reference, shm_name, index,dims,dtype):
    existing_shm = shared_memory.SharedMemory(name = shm_name)
    data_arr=np.ndarray(dims,dtype = dtype, buffer = existing_shm.buf)
    np_src = data_arr[:,:,index]
    output = np.array([])
    ref_data = reference.data                               ####OPTIMIZATION POSSIBLE HERE####
    try:
        output = astroalign.register(np_src, ref_data)[0]
    except KeyError:
        np_src = np_src.byteswap().newbyteorder()
        output = astroalign.register(np_src, ref_data)[0]
    data_arr[:,:,index] = output[:,:]
    existing_shm.close()


def align(hduls, read_ext="SCI", write_ext="ALGN", ref=None):
    """
    Aligns the source astronomical image(s) to the reference astronomical image
    \b
    :param hduls: list of fitsfiles
    :param name: header containing image data.
    :param reference: index of refence image. Should be file with greatest SNR
    :return: list of FITS files with <read_ext> HDU aligned
    """

    faulthandler.enable()
    snr_arr = np.zeros(len(hduls), dtype = np.float64)
    shm_snr = shared_memory.SharedMemory(create = True, size = snr_arr.nbytes)
    arr = np.ndarray(snr_arr.shape, dtype = snr_arr.dtype, buffer=shm_snr.buf)
    hduls = copy.deepcopy(hduls)
    snr(shm_snr.name, hduls, read_ext)
    snr_arr = np.zeros(len(hduls), dtype = np.float64)
    for i in range(len(hduls)):
        snr_arr[i] = arr[i]
    # No reference index given. we establish reference based on best signal to noise ratio
    if ref is None:
        try:
            reference = hduls[0][read_ext]  # 0th index reference is used by default
        except KeyError:
            reference = hduls[0]["PRIMARY"]
        ref_snr = arr[0]

        for index in range(len(hduls)): # loops though hdul lists finding the hdul with greatest snr
            if snr_arr[index] > ref_snr:  # compares SNR value of current hdul to refyee
                ref_snr = snr_arr[index]
                try:
                    reference = hduls[index][read_ext]
                except KeyError:
                    reference = hduls[index]["PRIMARY"]
    else:  # ref index is provided
        try:
            reference = hduls[ref][read_ext]
        except KeyError:
            reference = hduls[ref]["PRIMARY"]

    try:
        ref_data = reference.data
    except AttributeError:
        print("The reference file have doesn't have Attribute: Data")

    processes = []
    imagesize = ref_data.shape
    length = len(hduls)
    dims = (imagesize[0],imagesize[1],length) 
    data_array = np.ndarray(dims,dtype = np.float32)
    begin = datetime.datetime.now()
    for index in range(len(hduls)):
        try:
            data_array[:,:,index] = hduls[index][read_ext].data
        except KeyError:
            data_array[:,:,index] = hduls[index]["PRIMARY"].data
    shm_hdul_data = shared_memory.SharedMemory(create = True,size = data_array.nbytes)
    arr = np.ndarray(data_array.shape,dtype = data_array.dtype,buffer = shm_hdul_data.buf)
    arr[:,:,:] = data_array[:,:,:]
    for index in range(len(hduls)):
        p = mp.Process(target = multiprocessed_align, args = (reference, shm_hdul_data.name, index, dims, data_array.dtype))
        processes.append(p)
        p.start()
    for i in processes:
        i.join()
    image = np.zeros(imagesize,dtype = np.float32)
    data_array[:,:,:] = arr[:,:,:]
    for i in range(len(hduls)):
        image = data_array[:,:,i]
        try:
            hduls[i][read_ext].data = image
            idx = hduls[i].index_of(read_ext)
        except KeyError:
            hduls[i]["PRIMARY"].data = image
            idx = hduls[i].index_of("PRIMARY")
        hduls[i][idx].header['EXTNAME'] = (write_ext)
        hduls[i][idx].header = reference.header  
    print('\n' + 'Time for aligning {} images was {} seconds meaning {} seconds per image'.format(len(hduls), datetime.datetime.now()-begin, (datetime.datetime.now()-begin)/len(hduls)))
    shm_hdul_data.close()
    shm_snr.close()
    shm_snr.unlink()
    shm_hdul_data.unlink()
    return hduls


@cli.cli.command("align")
@click.option("-r", "--read_ext", default="SCI", help="The HDU to be aligned.")
@click.option("-w", "--write_ext", default="ALGN", help="The extension name that the resulting HDUL gets written to. Default""is `ALGN`")
@cli.operator
# TODO: use CAT sources if they exist

# align function wrapper
def align_cmd(hduls, read_ext="SCI", write_ext="ALGN", ref=None):
    """
    Aligns the source astronomical image(s) to the reference astronomical image
    \b
    :param hduls: list of fitsfiles
    :return: list of fistfiles with <read_ext> HDU aligned
    """
    return align([hduls for hduls in hduls], read_ext, write_ext, ref)