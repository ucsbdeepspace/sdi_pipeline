import tkinter as tk
from tkinter import filedialog
from astropy.io import fits
from glob import glob
import sys
import os

sdi_path = os.environ['HOME'] + "/sdi_pipeline/sdi"
sys.path.append(sdi_path) #Find a better way @Varun
from secid import secid

global sci
#Hides the Tk() window
root = tk.Tk()
root.withdraw()
 
#Opens the file explorer for directory selection
directory = filedialog.askdirectory(initialdir = "/home/alex/smalldata/science_data", title = "Select a Directory")

if directory == ():  #If the user does not select a file, close the script 
    print("File not selected. Closing script...")
    quit()

fz_glob = directory + "/*.fz"
fits_glob = directory + "/*.fits"

filenames = glob(fz_glob) + glob(fits_glob)
if filenames == []:
    raise ValueError('No fitz files in the selected directory')

for i,filename in enumerate(filenames):
    image = fits.open(filename)
    secid_gen = secid(image, "SCI")
    #print(secid_gen)
    hdul = next(secid_gen)
    print(f"{filename} is in section {hdul['SCI'].header['SECID']}")
    #print(output)

