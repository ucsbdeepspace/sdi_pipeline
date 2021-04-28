"""
light_curves finds the intensity of a list of sources over time
History:
    Created on 2019-09-06
e       Varun Iyer <varun_iyer@ucsb.edu>
"""
from ast import literal_eval
from pyds9 import DS9
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from . import align
from .sources import db


def image(hdul, cat, num=20):
    d = DS9()
    d.set_pyfits(hdul)
    for x, y in zip(cat['x'], cat['y']):
        d.set('regions command {{circle {} {} 40 #text=""}}'.format(x, y))


def curves(record_s):
    """
    Generate a plot of light curves for each source
    Arguments:
        sources -- list of catalog file HDUs
        images -- list of science images
    Keyword Arguments:
        num=100 -- how many of the curves to graph
        detected=100 -- list of sources to graph in a different color as detected
        fname="" -- path to save the image; if not "", the image is saved to the
            specified path
        show=False -- whether to display the plot
    """
    records = record_s
    if not isinstance(record_s, list):
        records = [record_s]

    for record in records:
        sources = sorted(record.sources, key=lambda x: x.image.time)
        times = [s.image.time for s in sources]
        fluxes = [literal_eval(s.data)[7] for s in sources]
        plt.plot(times, fluxes)
    plt.gcf().autofmt_xdate()
    plt.show()
