from setuptools import setup, find_packages

setup(
    name="sdi-cli",
    version="0.99",
    py_modules=["sdi"],
    # packages=find_packages(include=["openfits"]),
    include_package_data=True,
    setup_requires = ['numpy', "setuptools_rust"]
    install_requires=["click", "astropy", "photutils", "ois", "astroalign", "astroquery", "sklearn"],
    entry_points="""
        [console_scripts]
        sdi=sdi._cli:cli
    """,
)
