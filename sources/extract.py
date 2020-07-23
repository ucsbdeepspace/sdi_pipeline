"""
Extracts variable sources
History:
    Created by 2019-09-05
        Varun Iyer <varun_iyer@ucsb.edu>
"""
import os
import sep
from skimage.feature import blob_log
from ..common import to_np


def extract(residual_s, thresh=None, method="sep"):
    """
    Uses sep to find sources on a residual image(s)
    Arguments:
        residuals -- image of residuals from hotpants or a list of images
        method    -- algorithm to use, sep or skimage
    Returns:
        A list of Source objects representing the location and various metrics
            of detected variable sources
    """
    residuals = []
    if isinstance(residuals, list):
        residuals = residual_s
    else:
        residuals.append(residual_s)
    sources = []
    for r in residuals:
        r_np = to_np(r)
        if thresh is None:
            # from astroalign’s settings
            if method=="sep":
                print("sep")
                r_np = r_np.byteswap().newbyteorder()
                bkg = sep.Background(r_np)
                sources.append(sep.extract(r_np - bkg.back(), bkg.globalrms * 3.0, segmentation_map=False)) # compare overall rms
            if method=="skimage":
                print("skimage")
                sources.append(blob_log(r_np, min_sigma=1, max_sigma=20, num_sigma=10, threshold=0.2))
    return sources if isinstance(residuals, list) else sources[0]
