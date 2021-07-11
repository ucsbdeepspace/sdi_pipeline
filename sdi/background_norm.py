# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 19:01:30 2021

@author: prera

IN ORDER TO RUN THIS DOWNLOAD THE FOLLOWING files:
collate.py
"""

from glob import glob
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

#Requires collate - pathing currently incorrect
from collate import hdultocluster,cluster_search

import scipy.optimize as scp
from astropy.time import Time

#Organising data
images_path = glob(os.path.join('section_182','*.fz'),recursive=True)
images = [fits.open(p) for p in images_path]

times = [im['SCI'].header['DATE-OBS'] for im in images]
t = Time(times)
times = t.mjd

path = glob(os.path.join('sditest_sec182','*.fits'),recursive=True)
hduls = [fits.open(p) for p in path]

#xrts is a list of 99 catalogs files containing info about variables stars
xrts = [h['XRT'].data for h in hduls]

#The catalog file that identifies similar sources:
cat_files = glob(os.path.join('section_182_collated','*.fits'),recursive=True)
cats = [fits.open(f) for f in cat_files]
#objects is a set 517 catalog files, each of which refer to a single source.
#Each file then has 99 sets of data, that refers to the information about the star across the 99 images.
objects = hdultocluster(cats)


#find the variable stars based on one image
indices_var = []
for i in range(0,len(xrts[0])):
    indices_var.append(cluster_search(objects,xrts[0]['x'][i],xrts[0]['y'][0]))

indices_var = np.sort(np.array(list(dict.fromkeys(indices_var))))
indices_bg = [i for i in range(len(objects))]
for i in indices_var:
    indices_bg.remove(i)

#Finding the bg_stars:
bg_stars = [objects[i] for i in indices_bg]
#Find the ones in every image:
intermediate = []
for star in bg_stars:
    if star['flux'].all() != 0:
        intermediate.append(star)
bg_stars = intermediate

"""Now we are done collecting the data we need for analysis"""

#Plot background stars' flux
#for i in intermediate:
#    plt.hist(i['flux'])
#    plt.show()
plt.scatter(bg_stars[20]['flux'],bg_stars[30]['flux'])
#plt.scatter(bg_stars[30]['flux'],times)
plt.show()

#Based on this plot, we see that flux is distributed around a gaussian. Not too many outliers, so we won't do outlier analysis.

"fluxes are all the same, so y = x + gaussian error"

#Outlining a gaussian model where y will be our measurements fluxes, and x will be the exact flux we want to find
#The following function comes from the Homework 3 solutions by Professor Brandt, modified by me

def loglike_flux(p, table):

    norm_flux = p
    flux = table['flux']
    var_flux = table['fluxerr']**2
    
    like_member = np.exp(-(flux - norm_flux)**2/var_flux)
    like_member /= np.sqrt(2*np.pi*var_flux)
    
    return -np.sum(np.log(like_member))

tables = [bg_stars[i] for i in range(len(bg_stars))]
means = [np.mean(bg_stars[i]['flux']) for i in range(len(bg_stars))]
flux_normal = [scp.minimize(loglike_flux, means[i], (tables[i]), method="Powell").x for i in range(len(tables))]

#Calculating the confidence interval: Errors = 1sigma/sqrt(N)
errors = [np.std(bg_stars[i]['flux'])/np.sqrt(len(bg_stars[i])) for i in range(len(bg_stars))]

"""Now to normalize the variable stars w.r.t. the background stars we found
Ideally, we would know what variable star we are looking at, 
look up its color in a reference table and choose the background stars that have similar colours,
however in our case we do not know the colour of the variable source, so I will not be dealing with colour.
If we know what variable star we are looking at, we would also be able to form a model for its variability and be able to normalize from there
However, with no further information about the star we would not be able to do any sort of normalization of variable stars."""

#finding the variable stars
var_stars = [objects[i] for i in indices_var]
#Find the ones in every image:
intermediate_var = []
for star in var_stars:
    if star['flux'].all() != 0:
        intermediate_var.append(star)
var_stars = intermediate_var
"""If we had a model for the star, we would be able to build a likelihood function and linearize and normalize the function
similar to above, or use gaussian process region to buid a likelihood function."""

