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
        self.assertEqual(len(TestSubtract.output), 2042, "Did not yield 2042 subtractions")

    def test_output(self):
        # sorts the known true output to be congruent with the current output
        self.path_len = len(os.path.join(os.path.dirname(__file__), "fixtures/comparitiveData/subtractData"))
        self.paths = glob.glob("{}/*.fits*".format(os.path.join(os.path.dirname(__file__), "fixtures/comparitiveData/subtractData")))
        self.paths = sorted(self.paths, key = lambda item: int(item[self.path_len+1:len(item)-5]))

        # compares the knwon true output to current output
        self.true_output = []
        for p, o in zip(self.paths, TestSubtract.output):
            self.true_output = fits.open(p)
            compare = fits.FITSDiff(self.true_output, o)
            self.assertEqual(compare.identical, True , compare.report(fileobj = None))
            self.true_output.close()

    def test_click(self):
        runner = CliRunner()
        working = True
        try:
            runner.invoke(sdi._subtract_cmd)
        except:
            working = False
        self.assertEqual(working, True, "Click for subtract command not working")


if __name__ == "__main__":
    unittest.main()

