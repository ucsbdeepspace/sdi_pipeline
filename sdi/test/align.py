import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import fitsio
import sdi
import numpy as np
from skimage.util.dtype import img_as_float

class TestAlign(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/smallData")
    ## Defining these variables within setUpClass functions assigns wrong dtype to output
    read = [s for s in sdi.read(path)]
    output = sdi.align(read)

    @classmethod
    def tearDownClass(cls):
        for hdul in cls.read:
            hdul.close()

    def test_type(self):
        for o in TestAlign.output:
            self.assertIsInstance(o, fits.HDUList, "Did not align type fits.HDUList")

    def test_length(self):
        n = 0
        try:
            while True:
                next(TestAlign.output)
                n += 1
        except StopIteration: pass
        self.assertEqual(n, 2, "Did not align 2 fits files")


if __name__ == "__main__":
    unittest.main()
