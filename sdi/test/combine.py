"""
Author: Benjamin Fogiel - bfogiel@ucsb.edu
Tests the combine function
"""

import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi


class TestCombine(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/science")

    @classmethod
    def setUpClass(cls):
        cls.read = [s for s in sdi.read(cls.path)]
        cls.output = sdi.combine(cls.read)
        true_path = os.path.join(
            os.path.dirname(__file__), "fixtures/comparativeData/combineData/0.fits"
        )
        cls.true_output = [fits.open(true_path)]

    @classmethod
    def tearDownClass(cls):
        for o, t in zip(cls.read, cls.true_output):
            o.close()
            t.close()

    def test_length(self):
        self.assertEqual(
            len(TestCombine.output), 1, "Combine did not return one PrimaryHDU"
        )

    def test_type(self):
        self.assertIsInstance(
            TestCombine.output[0], fits.PrimaryHDU, "Output is not of type PrimaryHDU"
        )

    def test_output(self):
        self.compare = fits.FITSDiff(
            fits.HDUList(TestCombine.output[0]),
            fits.HDUList(TestCombine.true_output[0]),
        )
        self.assertEqual(
            self.compare.identical, True, self.compare.report(fileobj=None)
        )

    def test_click(self):
        runner = CliRunner()
        result = runner.invoke(sdi._combine_cmd, ["-n", "name"])
        assert result.exit_code == 0


if __name__ == "__main__":
    unittest.main()
