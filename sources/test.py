from astropy.io import fits
from glob import glob
filenames = glob("/home/varun/smalldata/*")
hduls = [fits.open(f) for f in filenames]
sci = [h['SCI'] for h in hduls]
cat = [h['CAT'] for h in hduls]
