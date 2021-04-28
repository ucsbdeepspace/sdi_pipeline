import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi

class TestRead(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/science")
    nopath = os.path.join(os.path.dirname(__file__), "fixtures/gobbledygook")

    def setUp(self):
        self.output = [s for s in sdi.read(TestRead.path)]
        self.noutput = [s for s in sdi.read(TestRead.nopath)]

    def tearDown(self):
        for hdul in self.output:
            hdul.close()

    def test_length(self):
        self.assertEqual(len(self.output), 10, "Did not read ten HDULs in.")
        self.assertEqual(len(self.noutput), 0, "Did not read zero HDULs in from " \
                                          "empty directory.")

    def test_type(self):
        for o in self.output:
            self.assertIsInstance(o, fits.HDUList, "Did not read type fits.HDUList")

    def test_tk(self):
        runner = CliRunner()
        result = runner.invoke(sdi._read_cmd)
        # FIXME figure out how to do click right

if __name__ == "__main__":
    unittest.main()
