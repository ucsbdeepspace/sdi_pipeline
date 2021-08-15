import tkinter as tk
from tkinter import filedialog
from astropy.io import fits
import sys
import os

sdi_path = os.environ['HOME'] + "/sdi_pipeline/sdi"
sys.path.append(sdi_path) #Find a better way @Varun
from secid import secid

#Hides the Tk() window
root = tk.Tk()
root.withdraw()
 
#Opens the file explorer for fits file selection
filenames = filedialog.askopenfilenames(initialdir = "/home/alex/smalldata/science_data", title = "Select a File", filetypes = (("Compressed Fitz files", "*.fz*"), ("Fitz files", "*.fits"), ("all files", "*.*")))

if filenames == ():  #If the user does not select a file, close the script 
    print("File not selected. Closing script...")
    quit()

hduls = [fits.open(f) for f in filenames]
secid_gen = secid(hduls, "SCI")

for i,hdul in enumerate(secid_gen):
#    print(hdul['SCI'].header)
    print(f"{filenames[i]} is in section {hdul['SCI'].header['SECID']}")
    #print(output)
