import sys
import click
from . import _cli as cli
from astropy.io import fits
import glob
from tkinter.filedialog import askdirectory

def read(directory):
    """
    Takes a directory containing fits files and returns them as a list
    """
    paths = glob.glob("{}/*.fits*".format(directory))
    hduls = [fits.open(p) for p in paths]
    return hduls

def write(hduls, directory, format_):
    """
    Writes all given hduls into the directory specified
    """
    import os
    for i, h in enumerate(hduls):
        path = os.path.join(directory, format_.format(number=i))
        click.echo(f"writing hdul to {path}")
        h.writeto(path)
    return hduls

@cli.cli.command("write")
@click.option('-d', '--directory', type=str, help="Specify path to directory to save fitsfiles.", default="./")
@click.option('-f', '--format', "format_", type=str, help="Specify string format for filename.", default="{number}.fits")
@cli.operator
## write function wrapper
def write_cmd(hduls, directory, format_):
    """
    Writes all given hduls into the directory specified
    """
    return write(hduls, directory, format_)

@cli.cli.command("read")
@click.option('-d', '--directory', type=str, help="Specify path to directory of fitsfiles.", required=False)
@cli.generator
## read function wrapper
def read_cmd(directory):
    """
    Takes a directory containing fits files and returns them as a list
    """
    # add try (evaluate if the try is efficient)
    if directory is None:
        try:
            directory = askdirectory()
        except:
            click.echo("Visual file dialog does not exist, please use option -d and specify path to directory to read fitsfiles.", err=True)
            sys.exit()
    hduls = read(directory)
    if hduls:
        return hduls
    else:
        sys.exit(f"Could not open fitsfiles from directory {directory}")
