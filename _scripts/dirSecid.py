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
directory = filedialog.askdirectory(title = "Select a Directory")

if directory == ():  #If the user does not select a directory, close the script 
    print("Directory not selected. Closing script...")
    quit()

fz_glob = directory + "/*.fz"
fits_glob = directory + "/*.fits"

filenames = glob(fz_glob) + glob(fits_glob)
if filenames == []:
    raise ValueError('No fitz files in the selected directory')

hduls = [fits.open(f) for f in filenames]
for i,hdul in enumerate(secid(hduls,"SCI")):
    print(f"{filenames[i]} is in section {hdul['SCI'].header['SECID']}")
    #print(output)

