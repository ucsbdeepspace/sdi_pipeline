from setuptools import setup, find_packages
import subprocess
import sys
import os

'''
try:
    import numpy
except ModuleNotFoundError:
    import sys
    sys.exit("numpy not found, tripp requires numpy for installation.\n Please try '$pip3 install numpy'.")

try:
    import setuptools_rust
except ModuleNotFoundError:
    import sys
    sys.exit("setuptools_rust not found, tripp requires setuptools_rust for installation.\n Please try '$pip3 install setuptools_rust'.")
'''

subprocess.run(sys.executable + " -m pip install --upgrade pip setuptools==56.0.0 setuptools_rust numpy", shell=True)
os.system('cp ./tripp/tripp/ois.py ./lib/python3.9/site-packages/')

setup(
    name="tripp-cli",
    version="0.99",
    py_modules=["tripp"],
    # packages=find_packages(include=["openfits"]),
    include_package_data=True,
    install_requires=["click", "astropy", "matplotlib", "regions", "photutils", "ois", "numba", "cupy-cuda11x", "astroalign", "astroquery", "scikit-learn"],
    entry_points="""
        [console_scripts]
        tripp=tripp._cli:cli
    """,
)
