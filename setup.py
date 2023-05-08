from setuptools import setup, find_packages
import subprocess
import sys

'''
try:
    import numpy
except ModuleNotFoundError:
    import sys
    sys.exit("numpy not found, sdi requires numpy for installation.\n Please try '$pip3 install numpy'.")

try:
    import setuptools_rust
except ModuleNotFoundError:
    import sys
    sys.exit("setuptools_rust not found, sdi requires setuptools_rust for installation.\n Please try '$pip3 install setuptools_rust'.")
'''

subprocess.run(sys.executable + " -m pip install --upgrade pip setuptools==56.0.0 setuptools_rust numpy", shell=True)

setup(
    name="sdi-cli",
    version="0.99",
    py_modules=["sdi"],
    # packages=find_packages(include=["openfits"]),
    include_package_data=True,
    install_requires=["cupy-cuda11x", "click", "astropy", "photutils", "astroalign", "astroquery", "scikit-learn","numpy==1.23.5","numba==0.53.1"],
    entry_points="""
        [console_scripts]
        sdi=sdi._cli:cli
    """,
)
