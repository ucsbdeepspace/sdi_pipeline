from astropy.io import fits
from sys import argv
from glob import glob
from . import align
from .combine import combine
from .subtract import subtract
from .sources import extract
import pickle

file_list = argv[1:]
if len(file_list) < 2:
    print("ERROR: Not enough filenames specified. Try\n\t python3 -m sdi_pipeline /path/to/fitsfiles/*.fz\n or something similar.")
    sys.exit(1)
science_images = [fits.open(f)["SCI"] for f in file_list]
aligned = align.image(science_images)
print("Finished alignment")
template = combine(aligned)
print("Finished combine")
residuals = subtract(aligned, template)
print("Finished subtract")

im_sources = extract([r[0] for r in residuals])

pickle.dump(im_sources, open("transient_candidates.pkl","wb"))

print("Finished extract")
