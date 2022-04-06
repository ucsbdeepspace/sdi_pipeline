"""
Unit testing code taken from Benjamin Fogiel
Tests collate command
"""

import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi
import glob


class TestCollate(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/extractedsims")

    @classmethod
    def setUpClass(cls):
        cls.read = [s for s in sdi.read(cls.path)]
        cls.output = list(sdi.collate(cls.read))
        cls.paths = glob.glob("fixtures/comparativeData/collateData/*")
        cls.true_output = [fits.open(path) for path in cls.paths]

    @classmethod
    def tearDownClass(cls):
        for r, t in zip(cls.read, cls.true_output):
            r.close()
            t.close()

    def test_type(self):
        for o in TestCollate.output:
            self.assertIsInstance(o, fits.HDUList, "Output is not type fits.HDUList")

    def test_length(self):
        self.assertEqual(
            len(TestCollate.output),
            10,
            "Did not collate correct number of HDULs from fixtures/extractedsims",
        )

    def test_output(self):
        for t, o in zip(TestCollate.true_output, TestCollate.output):
            compare = fits.FITSDiff(fits.HDUList(t), fits.HDUList(o))
            self.assertEqual(compare.identical, True, compare.report(fileobj=None))

    def test_click(self):
        runner = CliRunner()
        result = runner.invoke(sdi._collate_cmd, ["-c", "xy"])
        assert result.exit_code == 0


if __name__ == "__main__":
    unittest.main()
