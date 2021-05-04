import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi

class TestExtract(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/residuals")

    @classmethod
    def setUpClass(cls):
        cls.read = [s for s in sdi.read(cls.path)]
        cls.output = sdi.extract(cls.read)

    @classmethod
    def tearDownClass(cls):
        for hdul in cls.read:
            hdul.close()

    def test_length(self):
        self.assertEqual(len(TestExtract.output), 10, "Extract did not return 10 PrimaryHDU files")

    def test_length(self):
        n = 0
        try:
            while True:
                next(TestExtract.output)
                n += 1
        except StopIteration: pass
        self.assertEqual(n, 10, "Did not extract 10 fits files")

if __name__ == "__main__":
    unittest.main()

