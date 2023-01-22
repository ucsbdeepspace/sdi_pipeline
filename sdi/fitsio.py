"""
fitsio -- this module contains the read and write function
"""

import glob
from tkinter.filedialog import askdirectory
import sys
import click
from astropy.io import fits
import multiprocessing as mp
from . import _cli as cli
import os
import time

def mp_write(directory, hdul, i,format_):
    path = os.path.join(directory, format_.format(number=i))
    hdul.writeto(path)

def read(directory):
    """
    Takes a directory containing fits files and returns them as a list
    """
    start_time = time.perf_counter()
    paths = glob.glob("{}/*.fits*".format(directory))
    try:
        # if file name is numeric, paths will be sorted in ascending order
        paths = sorted(paths, key = lambda item: int(item[len(directory):len(item)-5]))
    except:
        pass
    hduls = [fits.open(p) for p in paths]
    stop_time = time.perf_counter()
    module_time = stop_time - start_time
    print(module_time)
    return hduls

def write(hduls, directory, format_):
    """
    Writes all given hduls into the directory specified
    """
    start_time = time.perf_counter()
    #check if directory exists
    is_dir = os.path.isdir(directory)

    if is_dir == False:
        os.mkdir(directory)
        print('directory was not present, now made at ' + directory)

    jobs=[]
    for i, h in enumerate(hduls):
        p = mp.Process(target = mp_write, args=(directory,h,i,format_))
        p.start()
        jobs.append(p)
#        path = os.path.join(directory, format_.format(number=i))
#        click.echo(f"writing hdul to {path}")
#        h.writeto(path)
    for j in jobs:
        j.join()
    stop_time = time.perf_counter()
    module_time = stop_time - start_time
    print(module_time)
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
