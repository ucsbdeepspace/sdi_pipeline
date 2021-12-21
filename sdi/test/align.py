"""
Author: Benjamin Fogiel - bfogiel@ucsb.edu
Tests the align function
"""

import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi
import numpy as np
from skimage.util.dtype import img_as_float
import glob

class TestAlign(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # uses two fits files as align param to decrease runtime
        cls.paths = glob.glob(os.path.join(os.path.dirname(__file__), "fixtures/science/tfn0m414-kb99-20180831-029[0-1]-e91.fits.fz"))
        cls.read = [fits.open(p) for p in cls.paths]
        cls.output = list(sdi.align(cls.read))

        # sorts the known true output to be congruent with the current output
        cls.path_len = len(os.path.join(os.path.dirname(__file__), "fixtures/comparativeData/alignData"))
        cls.paths_true = glob.glob("{}/*.fits*".format(os.path.join(os.path.dirname(__file__), "fixtures/comparativeData/alignData")))
        cls.paths_true = sorted(cls.paths_true, key = lambda item: int(item[cls.path_len+1:len(item)-5]))
        cls.true_output = [fits.open(p) for p in cls.paths_true]

    @classmethod
    def tearDownClass(cls):
        for r, t in zip(cls.read, cls.true_output):
            r.close()
            t.close()

    def test_type(self):
        for o in TestAlign.output:
            self.assertIsInstance(o, fits.HDUList, "Did not align type fits.HDUList")

    def test_length(self):
        self.assertEqual(len(TestAlign.output), 2, "Did not align the correct number of fits")

    def test_output(self):
        for t, o in zip(TestAlign.true_output, TestAlign.output):
            ## converts the true_output SNR value to type float32 to be congruent with output
            t[0].header['SNR'] = np.float32(t[0].header['SNR'])
            compare = fits.FITSDiff(t, o)
            self.assertEqual(compare.identical, True , compare.report(fileobj = None))

    def test_click(self):
        runner = CliRunner()
        result = runner.invoke(sdi._align_cmd, ['-n', 'name'])
        assert result.exit_code == 0


if __name__ == "__main__":
    unittest.main()
