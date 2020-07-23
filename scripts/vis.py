from pyds9 import DS9
import numpy as np

#hduls should be a list of hdul(s), cat should be a string containing the path to a .npy file
#Adding circles via the image() function will apply circles to all images
def image(hduls, cat=None, color='green', size=40):
    d = DS9()
    d.set("scale mode zscale")
    print("asdf")

    if cat == None:
        for hdul in hduls:
            d.set("frame new")
            d.set_pyfits(hdul)
            d.set("zoom to fit")
    else:
        print(len(cat))
        print(len(hduls))
        for c, hdul in zip(cat, hduls):
            print("stuff!")
            d.set("frame new")
            d.set_pyfits(hdul)
            d.set("zoom to fit")
            # coords = np.load(cat)
            for coord in c:
                d.set('regions command {{circle {} {} {} #color={}}}'.format(coord[0], coord[1], size, color))
    d.set("frame first")
    d.set("frame delete")

#Use this function to add circles to currently selected image in DS9
def circle(cat, color='green', size=40):
    coords = np.load(cat)
    d = DS9()
    for coord in coords:
        d.set('regions command {{circle {} {} {} #color={}}}'.format(coord[0], coord[1], size, color))
