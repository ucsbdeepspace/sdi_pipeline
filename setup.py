from setuptools import setup, find_packages

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
