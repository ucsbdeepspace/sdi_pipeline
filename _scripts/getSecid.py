import tkinter as tk
from tkinter import filedialog
from astropy.io import fits
import sys
import os

sdi_path = os.environ['HOME'] + "/sdi_pipeline/sdi"
sys.path.append(sdi_path) #Find a better way @Varun
from secid import secid

global sci
#Hides the Tk() window
root = tk.Tk()
root.withdraw()
 
#Opens the file explorer for fits file selection
path = filedialog.askopenfilename(initialdir = "/home/alex/smalldata/science_data", title = "Select a File", filetypes = (("Compressed Fitz files", "*.fz*"), ("Fitz files", "*.fits"), ("all files", "*.*")))

if path == ():  #If the user does not select a file, close the script 
    print("File not selected. Closing script...")
    quit()

image = fits.open(path)

#image = fits.open(path)['PRIMARY']

secid_gen = secid(image, "SCI")
#print(secid_gen)
output = next(secid_gen)
print(f"{path} is in section {output[0].header['SECID']}")
#print(output)

#for value in secid_gen:
#    print(value)
