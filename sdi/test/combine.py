import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi

class TestCombine(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/science")
#    nopath = os.path.join(os.path.dirname(__file__), "fixtures/gobbledygook")

    def setUp(self):
        self.read = [s for s in sdi.read(TestCombine.path)]
        self.output = sdi.combine(self.read)

        #        self.noutput = [s for s in sdi.read(TestRead.nopath)]

    def tearDown(self):
        for hdul in self.read:
            hdul.close()

    def test_length(self):
        self.assertEqual(len(self.output), 1, "Did not combine HDULs into one PrimaryHDU")
#        self.assertEqual(len(self.noutput), 0, "Did not read zero HDULs in from " \
#                                          "empty directory.")

    def test_type(self):
        self.assertIsInstance(self.output[0], fits.PrimaryHDU, "Output is not of type PrimaryHDU")

    def test_tk(self):
        runner = CliRunner()
        result = runner.invoke(sdi._read_cmd)
        # FIXME figure out how to do click right

if __name__ == "__main__":
    unittest.main()

