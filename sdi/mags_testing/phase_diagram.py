import numpy as np
from astropy.timeseries import LombScargle
from PyAstronomy.pyasl import foldAt
from astropy.time import Time
#from collate import hdultocluster,cluster_search
from glob import glob
import os
from astropy.io import fits
from astropy.time import Time
from scipy import signal
import matplotlib.pyplot as plt

''' Let's define a class that takes in input lightcurves
performs a lombscargle on them to find the frequencies
that we need, and phase folds on the interesting
frequencies'''
class lightcurve():
    # This function (constuctor) takes in user input
    # of the lightcurve time, flux and error.
    def __init__(self, time, flux, err, freq=None):
        self.time = time
        self.flux = flux
        self.err = err
#        self.N = len(time) # Number of points in the lightcurve
        # We will assign these values after the lombscargle
        self.period = None
        self.freq = freq
        self.power = None
        # We will assign these values after phase-folding
        self.phase = None
        self.folded_flux = None
    # Take the lombscargle and set the lightcurve period
    # to be the maximum of the lombscargle
    def lombscargle(self):
        if self.freq is None:
           self.freq, self.power = LombScargle(self.time, self.flux, self.err).autopower(minimum_frequency=1,maximum_frequency =3)
        else:
           self.power = LombScargle(self.time, self.flux, self.err).power(self.freq)
        # Set the lightcurve period to be the frequency corresponding to the maximum power\n",
        self.period = 1 / self.freq[np.argmax(self.power)]
    # Generate the phase folded lightcurve
    def phase_fold(self):
        if not self.period:
            self.lombscargle()
        # Generate phases corresponding to each point in time
        # We are phase-folding on the period with max power in lombscargle
        phase = foldAt(self.time, self.period, T0 = self.time[0])
        # Get the sorted indices of the phase
        index = np.argsort(phase)
        # Sort both the flux and phase
        self.phase = np.array(phase)[index]
        self.folded_flux = np.array(self.flux)[index]

data = np.genfromtxt('/home/pkotta/CSV/lc_GTAnd.csv',delimiter='\t')
time = []
magauto = []
magerr = []
for row in data:
    time.append(row[0])
    magauto.append(row[1])
    magerr.append(row[2])
time_PTF = np.array(time)
magauto = np.array(magauto)
magerr_PTF = np.array(magerr)
my_freq = np.linspace(1.5,1.8,3000)
#frequency is 1.727652065
lc_PTF = lightcurve(time_PTF, magauto,magerr_PTF,my_freq)
lc_PTF.lombscargle()
lc_PTF.phase_fold()

data_SDI = np.genfromtxt('var_bright_test.csv',delimiter=',')
t = []
mag = []
magerr_SDI = []
for row in data_SDI:
    t.append(row[2])
    mag.append(row[0])
    magerr_SDI.append(row[1])
t = np.array(t)
mag = np.array(mag)
magerr = np.array(magerr_SDI)
my_freq = np.linspace(1.5,1.8,3000)
#frequency is 1.727652065
lc_SDI = lightcurve(t, mag,magerr,my_freq)
lc_SDI.lombscargle()
lc_SDI.phase_fold()
"""Data from our images"""
"""
hduls = glob(os.path.join('sditest_output','*.fits'),recursive=True)
im_files = [fits.open(h) for h in hduls]
#Pipeline output hduls should be run through collate
cat_files = glob(os.path.join(''),recursive=True)
cats = [fits.open(f) for f in cat_files]
objects = hdultocluster(cats)

#We can manually locate the x and y coordinates of the star on the image,
x_var = 315
y_var = 253

rr_star = cluster_search(objects,x_var,y_var)
t = [im_files[i][0].header['DATE-OBS'] for i in range(len(im_files))]
times = Time(t)
t = times.mjd
#Transform time to mjd maybe
mag = rr_star['mag']
magerr = rr_star['magerr']
"""
#setting up the window function. The duration of each exposure is much shorter than the time over which the images are taken so I am approximating to delta functions
"""
def delta_funcs(x, xmin=None, xmax=None):
   # Return arrays for plotting delta functions
   # at locations x with heights h
    h = 1
    if xmin is None:
        xmin = min(x) - 1
    if xmax is None:
        xmax = max(x) + 1
    dx = 0.002 * (xmax - xmin)
    def vals():
        yield (xmin, 0)
        for xi, hi in sorted(np.broadcast(x, h)):
            yield from zip([xi - dx, xi, xi + dx], [0, hi, 0])
        yield (xmax, 0)
    return zip(*vals())

obs_window=delta_funcs(t,xmin = np.min(t),xmax = np.max(t))
window_func = []

for i in range(len(t)):
    window_func.append(next(obs_window))
   
plt.figure()
plt.xlabel('Time (MJD)')
plt.plot(window_func[0],window_func[1])
plt.show()

#Comparing the transform of our window function to the produced lombscargle periodogram
window_transform = LombScargle(window_func[0],window_func[1]).autopower(minimum_frequency=1,maximum_frequency =3)
plt.figure()
plt.plot(window_transform[0],window_transform[1])
plt.xlabel('frequency (Hz)')
plt.ylabel('Power')
plt.show()
"""
lc_SDI = lightcurve(t,mag,magerr,np.linspace(1.3,1.8,3000))
#print(len(lc_SDI.time),np.max(lc_SDI.flux),np.min(lc_SDI.flux))
lc_SDI.lombscargle()
lc_SDI.phase_fold()

#Plotting
print(lc_SDI.period)
print(lc_PTF.period)
plt.errorbar(lc_PTF.phase, lc_PTF.folded_flux,yerr = magerr_PTF,fmt = 'r*')
plt.errorbar(lc_SDI.phase, lc_SDI.folded_flux, yerr =None ,fmt = 'b*')
#plt.plot(lc_PTF.time,lc_PTF.flux,'o')
#plt.plot(lc_SDI.time,lc_SDI.flux,'o')
plt.show()
plt.plot(lc_SDI.freq,lc_SDI.power)
plt.show()
plt.plot(lc_PTF.freq,lc_PTF.power)
plt.show()
