import unittest
import click
from . import _cli as cli
from . import test

@cli.cli.command("test")
@click.option('-n', '--name', type=str, help="Specify the name of the function you'd like to test", default="test")
@cli.operator
def test_cmd(hduls, name):
    """
    Runs all tests through unittest.
    """
    if name == "align":
    	test_suite = unittest.TestLoader().loadTestsFromModule(test.align)
    elif name == "extract":
    	test_suite = unittest.TestLoader().loadTestsFromModule(test.extract)
    elif name == "combine":
    	test_suite = unittest.TestLoader().loadTestsFromModule(test.combine)
    elif name == "read":
    	test_suite = unittest.TestLoader().loadTestsFromModule(test.fitsio)
    elif name == "write":
    	test_suite = unittest.TestLoader().loadTestsFromModule(test.fitsio)
    elif name == "subtract":
    	test_suite = unittest.TestLoader().loadTestsFromModule(test.subtract)
    elif name == "test":
    	test_suite = unittest.TestLoader().loadTestsFromModule(test)
    else:
    	print("Specify a valid function you'd like to test:")
    	print("align, extract, combine, read, write, subtract")
    	print("Note: if no option is provided, it will run all tests")
    	return hduls;

    runner = unittest.TextTestRunner()
    runner.run(test_suite)
    return hduls
