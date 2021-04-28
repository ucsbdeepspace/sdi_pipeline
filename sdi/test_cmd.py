import unittest
import click
from . import _cli as cli
from . import test

@cli.cli.command("test")
@cli.operator
def test_cmd(hduls):
    """
    Runs all tests through unittest.
    """
    test_suite = unittest.TestLoader().loadTestsFromModule(test)
    runner = unittest.TextTestRunner()
    runner.run(test_suite)
    return hduls
