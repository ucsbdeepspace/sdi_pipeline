"""
Author: Benjamin Fogiel - bfogiel@ucsb.edu
Tests the subtract function
"""

import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import tripp
import glob

class TestSubtract(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # uses one fits file as subtract param to decrease runtime
        cls.path = glob.glob(os.path.join(os.path.dirname(__file__), "fixtures/science/tfn0m414-kb99-20180831-029[0-1]-e91.fits.fz"))
        cls.read = [fits.open(p) for p in cls.path]
        cls.output = list(tripp.subtract(cls.read))

        # sorts the known true output to be congruent with the current output
        cls.path_len = len(os.path.join(os.path.dirname(__file__), "fixtures/comparativeData/subtractData"))
        cls.paths = glob.glob("{}/*.fits*".format(os.path.join(os.path.dirname(__file__), "fixtures/comparativeData/subtractData")))
        cls.paths = sorted(cls.paths, key = lambda item: int(item[cls.path_len+1:len(item)-5]))
        cls.true_output = [fits.open(p) for p in cls.paths]

    @classmethod
    def tearDownClass(cls):
        for o, t in zip(cls.output, cls.true_output):
            o.close()
            t.close()

    def test_type(self):
        for o in TestSubtract.output:
            self.assertIsInstance(o, fits.HDUList, "Did not output type fits.HDUList")

    def test_length(self):
        self.assertEqual(len(TestSubtract.output), 2, "Did not output the correct number of HDULs")

    def test_output(self):
        for t, o in zip(TestSubtract.true_output, TestSubtract.output):
            compare = fits.FITSDiff(fits.HDUList(t), fits.HDUList(o))
            self.assertEqual(compare.identical, True, compare.report(fileobj = None))

    def test_click(self):
        runner = CliRunner()
        result = runner.invoke(tripp._subtract_cmd, ['-n', 'name'])
        assert result.exit_code == 0


if __name__ == "__main__":
    unittest.main()

