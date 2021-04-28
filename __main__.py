from astropy.io import fits
from sys import argv
from glob import glob
from . import align
from .combine import combine
from .subtract import subtract
from .sources import extract
import pickle

@click.command()
@click.option('-s', is_flag=True)
@click.argument('path')
@click.argument('save_path', nargs= -1) #nargs= -1 allows this argument to be optional but saves value as tuple
def run(s, path, save_path):
    save = False
    if path[-1] != "/":
        print("Please give the path to the files WITHOUT '*.fz' at the end.")
        print("Aborted!")
        exit(0)
    if s: 
        print("residuals will be saved")
        save_base = os.path.basename(glob("{}*.fz".format(path))[0])
        save = True
    for image in glob("{}*.fz".format(path)):
        file_list.append(image)
    science_images = [fits.open(f)["SCI"] for f in file_list]
    aligned = align.image(science_images)
    print("Finished alignment")
    template = combine(aligned)
    print("Finished combine")
    residuals = subtract(aligned, template)
    print("Finished subtract")
    if save:
        print("Saving Residuals")
        ri = [r[0] for r in residuals]
        save_path = os.path.expanduser(save_path)
        for count, p in enumerate(ri):
            fits.PrimaryHDU(p).writeto(save_file.format(spath = list(save_path)[0], index = count, original = save_base))
    im_sources = extract([res[0] for res in residuals])
    print("Finished extract")
    pickle.dump(im_sources, open("transient_candidates.pkl","wb"))

<<<<<<< HEAD
im_sources = extract([r[0] for r in residuals])
=======
>>>>>>> 0059029a093f2ca2cb95cd1e924160f64b2f2141

if __name__ == '__main__': 
	run()

