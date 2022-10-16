import click
import ois
#import zogy.py
from astropy.io.fits import CompImageHDU
from . import _cli as cli
from .combine import combine
import time

def subtract(hduls, name="ALGN", method: ("ois", "numpy", "sfft") ="ois"):
    """
    Returns differences of a set of images from a template image
    Arguments:
        hduls -- list of fits hdul's where the last image is the template 
        name -- name of the HDU to use image that the other images will be subtracted from
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    hduls = [h for h in hduls]
    outputs = []
    template_fits = combine(hduls, name)
    template =template_fits["PRIMARY"].data # this is a temporary HDUL containing 1 PrimaryHDU with combined data ## Align is the reference image here??
    print(" ")
    print("Beginning Subtraction")
    start = time.perf_counter()
    if method == "ois":
        print("Method = OIS")
        for i,hdu in enumerate(hduls):
            try:
                diff = ois.optimal_system(image=hdu[name].data, refimage=template, method='Bramich')[0] #Align is the input image here??
            except ValueError:
                diff = ois.optimal_system(image=hdu[name].data.byteswap().newbyteorder(), refimage=template.byteswap().newbyteorder(), method='Bramich')[0]
            hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i]['ALGN'].header, name = "SUB"))
            outputs.append(hdu)
    elif method == "numpy":
        print("Method = Numpy")
        for i,hdu in enumerate(hduls):
            diff = template - hdu[name].data
            hdu.insert(1,CompImageHDU(data = diff, header =  hduls[i]['ALGN'].header, name = "SUB"))
            outputs.append(hdu)
    elif method == "sfft":
        print("Method = SFFT")
       # Write temporary fits files to feed into SFFT subtract
	# TODO: Remove if SFFT subtract updates to accept data stacks instead of full fits files
		print(" ")
		print("Writing Temporary Fits Files")
		start2 = time.perf_counter()
		temp_image_fits_filenames = []
		for i, hdu in enumerate(hduls): #Write the science hduls into temporary fits files
			temp_filename = "temp_image{}.fits".format(i)
			hdu[name].writeto(temp_filename, overwrite = True)
			temp_image_fits_filenames.append(temp_filename)

		template_fits = combine(hduls, name) # Create the reference image
		template_fits.writeto("temp_ref.fits", overwrite = True) # Write the reference image into temporary fits file

		stop2 = time.perf_counter()
		print("Writing Complete")
		print("Time Elapsed: {} sec".format(round(stop2-start2, 4)))
		
		# Check file size of temporary fits files on disc
		size = 0
		for i, fits_name in enumerate(temp_image_fits_filenames): 
			size += os.path.getsize(fits_name)
		size += os.path.getsize("temp_ref.fits")
		print("Total size of temporary files is {} mb".format(size/1024**2))
		outputs = []

		#Subtract
		print(" ")
		print("Beginning SFFT Subtraction")
		start2 = time.perf_counter()
		print("Method = SFFT")
		#TODO Incorperate input masks for when we take real data
		for i, fits_name in enumerate(temp_image_fits_filenames):       
			sol, diff = Customized_Packet.CP(FITS_REF = "temp_ref.fits", FITS_SCI = fits_name, 
									FITS_mREF = 'temp_ref.fits', FITS_mSCI = fits_name,
									ForceConv = "REF", GKerHW = 4, BGPolyOrder = 2, KerPolyOrder = 2)    
			hdul = fits.open(fits_name)
			if np.isnan(np.sum(diff)) == True:
				raise ValueError("Residual contains NaN")
			else:
				hdul.insert(1,CompImageHDU(data = diff, header =  hdul['ALGN'].header, name = "SUB"))
			outputs.append(hdul)
		stop = time.perf_counter()
		print("Subtraction Complete")
		print("Time Elapsed: {} sec".format(round(stop-start, 4)))
		
	#Remove Temporary Fits files from disc    
		print("Removing Temporary Fits Files")
		for i, fits_name in enumerate(temp_image_fits_filenames):
			if os.path.exists(fits_name):
				os.remove(fits_name)
			else:
				print("{} does not exist".format(fits_name))
		if os.path.exists("temp_ref.fits"):
				os.remove("temp_ref.fits")
		else:
			print("temp_ref.fits does not exist")
		print("Removal Complete")
	else:
        raise ValueError(f"method {method} unknown!")
    stop2 = time.perf_counter()
    print("SFFT Subtraction Complete")
    print("Time Elapsed: {} sec".format(round(stop2-start2, 4)))
    return (hdul for hdul in outputs)

@cli.cli.command("subtract")
@click.option("-n", "--name", default="ALGN", help="The HDU to be aligned.")
@click.option("-m", "--method", default="ois", help="The subtraction method to use; ois or numpy (straight subtraction) or sfft (GPU accelerated).")
@cli.operator

## subtract function wrapper
def subtract_cmd(hduls, name="ALGN", method="ois"):
    """
    Returns differences of a set of images from a template image\n
    Arguments:\n
        hduls -- list of fits hdul's where the last image is the template\n
        name -- name of the HDU to use image that the other images will be subtracted from\n
        method -- method of subtraction. OIS is default (best/slow). Numpy is alternative (worse/quick)
    """
    return subtract(hduls, name, method)
