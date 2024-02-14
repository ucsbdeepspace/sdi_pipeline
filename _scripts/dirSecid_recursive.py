import tkinter as tk
from tkinter import filedialog
from astropy.io import fits
from glob import glob
import sys
import os
from tripp.secid import secid

#Hides the Tk() window
root = tk.Tk()
root.withdraw()
 
#Opens the file explorer for directory selection
directory = filedialog.askdirectory(title = "Select a Directory",initialdir="/mnt/annex/science_data")

if directory == ():  #If the user does not select a directory, close the script 
    print("Directory not selected. Closing script...")
    quit()

fz_glob = directory + "/**/*.fz"
fits_glob = directory + "/**/*.fits"
filenames = glob(fz_glob,recursive=True) + glob(fits_glob,recursive=True)

if filenames == []:
    raise ValueError('No fitz files in the selected directory')

for filename in filenames:
    fits_file = [fits.open(filename)]
    hdul = list(secid(fits_file))[0]
#    for i,hdul in enumerate([secid(fits.open(filename))]):
#        print(hdul) 
    print(f"{filename} is in section {hdul['SCI'].header['SECID']}")
