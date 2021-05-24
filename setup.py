from setuptools import setup, find_packages
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


setup(
    name="sdi-cli",
    version="0.99",
    py_modules=["sdi"],
    # packages=find_packages(include=["openfits"]),
    include_package_data=True,
    install_requires=["click", "astropy", "photutils", "ois", "pyds9", "astroalign", "astroquery"],
    entry_points="""
        [console_scripts]
        sdi=sdi._cli:cli
    """,
)
