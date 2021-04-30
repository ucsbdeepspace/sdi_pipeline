import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi
import glob

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

class TestWrite(unittest.TestCase):

    def setUp(self):
        self.h = sdi.read(os.path.join(os.path.dirname(__file__), "fixtures/science"))
        self.output = sdi.write(self.h, os.path.join(os.path.dirname(__file__), "fixtures/TestDir"), "{number}.fits")
        self.paths = glob.glob("{}/*.fits*".format(os.path.join(os.path.dirname(__file__), "fixtures/TestDir")))

    def tearDown(self):
        for hdul in self.h:
            hdul.close()
        for p in self.paths:
            os.remove(p)

    def test_type(self):
        for o in self.output:
            self.assertIsInstance(o, fits.HDUList, "Did not read type fits.HDUList")

    def test_dir(self):
        self.assertEqual(len(self.paths), 10, "Did not write ten HDULs to directory {self.paths}")
    

if __name__ == "__main__":
    unittest.main()
