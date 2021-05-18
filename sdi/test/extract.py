"""
Author: Benjamin Fogiel - bfogiel@ucsb.edu
Tests the extract function
"""

import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi
import glob

class TestExtract(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/residuals")

    @classmethod
    def setUpClass(cls):
        cls.read = [s for s in sdi.read(cls.path)]
        cls.output = list(sdi.extract(cls.read))

        # sorts the known true output to be congruent with the current output
        cls.path_len = len(os.path.join(os.path.dirname(__file__), "fixtures/comparitiveData/extractData"))
        cls.paths = glob.glob("{}/*.fits*".format(os.path.join(os.path.dirname(__file__), "fixtures/comparitiveData/extractData")))
        cls.paths = sorted(cls.paths, key = lambda item: int(item[cls.path_len+1:len(item)-5]))
        cls.true_output = [fits.open(p) for p in cls.paths]

    @classmethod
    def tearDownClass(cls):
        for r, t in zip(cls.read, cls.true_output):
            r.close()
            t.close()

    def test_type(self):
        for o in TestExtract.output:
            self.assertIsInstance(o, fits.HDUList, "Output is not of type fits.HDUList")

    def test_length(self):
        self.assertEqual(len(TestExtract.output), 10, "Did not extract 10 fits files")

    def test_output(self):
        for t, o in zip(TestExtract.true_output, TestExtract.true_output):
            compare = fits.FITSDiff(fits.HDUList(t), fits.HDUList(o))
            self.assertEqual(compare.identical, True, compare.report(fileobj = None))

    def test_click(self):
        runner = CliRunner()
        working = True
        try:
            runner.invoke(sdi._extract_cmd)
        except:
            working = False
        self.assertEqual(working, True, "Click for extract command not working")

if __name__ == "__main__":
    unittest.main()

