import os
import unittest
import click
from click.testing import CliRunner
from astropy.io import fits
import sdi

class TestSubtract(unittest.TestCase):

    path = os.path.join(os.path.dirname(__file__), "fixtures/smallData")
    
    @classmethod
    def setUpClass(cls):
        cls.read = [s for s in sdi.read(cls.path)]
        cls.output = sdi.subtract(cls.read)

    @classmethod
    def tearDownClass(cls):
        for hdul in cls.read:
            hdul.close()

    def test_type(self):
        for o in TestSubtract.output:
            self.assertIsInstance(o, fits.HDUList, "Did not subtract type fits.HDUList")
    
    def test_length(self):
        n = 0
        try:
            while True:
                next(TestSubtract.output)
                n += 1
        except StopIteration: pass
        self.assertEqual(n, 4084, "Did not subtract 4084 differences")

if __name__ == "__main__":
    unittest.main()

