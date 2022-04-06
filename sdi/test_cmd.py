import sys
import unittest
import click
from . import _cli as cli
from . import test


@cli.cli.command("test")
@click.option(
    "-n",
    "--name",
    type=str,
    help="Specify the name of the function you'd like to test",
    default="test",
)
@cli.operator
def test_cmd(hduls, name):
    """
    Runs all tests through unittest

    Avaible tests: fitsio (read and write), align, combine, extract, and subtract

    To run all tests, simply don't pass any parameters
    """
    if name == "test":
        test_suite = unittest.TestLoader().loadTestsFromModule(test)
    else:
        test_suite = unittest.TestLoader().loadTestsFromModule(getattr(test, name))

    runner = unittest.TextTestRunner()
    output = runner.run(test_suite).wasSuccessful()
    if not output:
        sys.exit(1)
    return hduls
