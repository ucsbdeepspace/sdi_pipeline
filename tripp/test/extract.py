"""
Author: Benjamin Fogiel - bfogiel@ucsb.edu
Tests the extract function
"""

import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import tripp
import glob

class TestExtract(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/residuals")

    @classmethod
    def setUpClass(cls):
        cls.read = [s for s in tripp.read(cls.path)]
        cls.output = list(tripp.extract(cls.read))

        # sorts the known true output to be congruent with the current output
        cls.path_len = len(os.path.join(os.path.dirname(__file__), "fixtures/comparativeData/extractData"))
        cls.paths = glob.glob("{}/*.fits*".format(os.path.join(os.path.dirname(__file__), "fixtures/comparativeData/extractData")))
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
        self.assertEqual(len(TestExtract.output), 10, "Did not extract the correct number of HDULs from fixtures/residuals")

    def test_output(self):
        for t, o in zip(TestExtract.true_output, TestExtract.output):
            compare = fits.FITSDiff(fits.HDUList(t), fits.HDUList(o))
            self.assertEqual(compare.identical, True, compare.report(fileobj = None))

    def test_click(self):
        runner = CliRunner()
        result = runner.invoke(tripp._extract_cmd, ['-t', 3.0, '-r', 0, '-w', 'XRT', 1])
        assert result.exit_code == 0

if __name__ == "__main__":
    unittest.main()

