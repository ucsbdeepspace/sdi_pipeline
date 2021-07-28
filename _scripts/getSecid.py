import tkinter as tk
from tkinter import filedialog
from astropy.io import fits
from sdi_pipeline.sdi.secid import secid
global sci
#Hides the Tk() window
root = tk.Tk()
root.withdraw()

#Opens the file explorer for fits file selection
path = filedialog.askopenfilename(initialdir = "/home/alex/smalldata/science_data", title = "Select a File", filetypes = (("Compressed Fitz files", "*.fz*"), ("Fitz files", "*.fits"), ("all files", "*.*")))
 
#print(path)

image = fits.open(path)

#image = fits.open(path)['PRIMARY']

secid_gen = secid(image, "SCI")
#print(secid_gen)
output = next(secid_gen)
print(output[0].header['SECID'])
#print(output)

#for value in secid_gen:
#    print(value)
