import numpy as np
from astropy.io import fits
from astropy import wcs
import pickle

template_image = fits.open("/home/williamwang/Downloads/lsc0m412-kb26-20190726-0290-e91.fits.fz")

w = wcs.WCS(template_image["SCI"].header)

source_list = pickle.load(open("Bramich_subtracted_sources.pkl",'rb'))

x = source_list[0][0]['x']
y = source_list[0][0]['y']
coordinates = np.stack((x,y),axis=-1)

xlist = []
ylist = []

print("The RA and DEC for (0,0) is", w.wcs_pix2world((np.array([[0,0]])),0))

print("The RA and DEC for (3054,0) is", w.wcs_pix2world((np.array([[3054,0]])),0))

print("The RA and DEC for (0,2042) is", w.wcs_pix2world((np.array([[0,2042]])),0))

print("The RA and DEC for (3054,2042) is", w.wcs_pix2world((np.array([[3054,2042]])),0))

for i in coordinates:
    pixarray = np.array([[i[0],i[1]]])
    radec = w.wcs_pix2world(pixarray,0)
    xlist.append(radec[0][0])
    ylist.append(radec[0][1])

xarray = np.array(xlist)
yarray = np.array(ylist)

radecs = np.array([xarray, yarray])

radecs = radecs.T

with open('source_radecs_section25.txt', 'w') as source_radecs:
    
    np.savetxt(source_radecs, radecs, delimiter=",", fmt=["%.8f", "%.8f"])
