"""
Author: Benjamin Fogiel - bfogiel@ucsb.edu
Tests the subtract function
"""

import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi
import glob

class TestSubtract(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # uses one fits file as subtract param to decrease runtime
        cls.path = glob.glob(os.path.join(os.path.dirname(__file__), "fixtures/science/tfn0m414-kb99-20180831-0290-e91.fits.fz"))
        cls.read = [fits.open(p) for p in cls.path]
        cls.output = list(sdi.subtract(cls.read))

    @classmethod
    def tearDownClass(cls):
        for o in cls.output:
            o.close()

    def test_type(self):
        for o in TestSubtract.output:
            self.assertIsInstance(o, fits.HDUList, "Did not output type fits.HDUList")

    def test_length(self):
        self.assertEqual(len(TestSubtract.output), 1, "Did not output the correct number of HDULs")

    def test_output(self):
        self.true_output = os.path.join(os.path.dirname(__file__), "fixtures/comparativeData/subtractData/0.fits")
        for o in TestSubtract.output:
            self.compare = fits.FITSDiff(fits.HDUList(o), self.true_output)
            self.assertEqual(self.compare.identical, True, self.compare.report(fileobj = None))


    def test_click(self):
        runner = CliRunner()
        result = runner.invoke(sdi._subtract_cmd, ['-n', 'name'])
        assert result.exit_code == 0


if __name__ == "__main__":
    unittest.main()

