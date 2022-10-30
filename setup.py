import os
import sys
from setuptools import setup, find_packages

DESCRIPTION = "Pipeline to Take a Stack of Images and Extract Residuals"
NAME = "sdi-pipeline"
AUTHOR = "Kyle Lam"
AUTHOR_EMAIL = "kylelam@ucsb.edu"
MAINTAINER = "Kyle Lama"
MAINTAINER_EMAIL = "kylelam@ucsb.edu"
DOWNLOAD_URL = "https://github.com/ucsbdeepspace/sdi_pipeline.git"

LICENSE = 'MIT Licence'
VERSION = '0.99'

install_reqs = ["numpy",
				"setuptools_rust",
				"click", 
				"astropy", 
				"photutils", 
				"ois", 
				"astroalign", 
				"astroquery", 
				"sklearn",
				]
setup_reqs = ["numpy", 
				"setuptools_rust"]

setup(name = NAME,
      version = VERSION,
      description = DESCRIPTION,
	  setup_requires = setup_reqs,
      install_requires = install_reqs,
      author = AUTHOR,
      author_email = AUTHOR_EMAIL,
      maintainer = MAINTAINER,
      maintainer_email = MAINTAINER_EMAIL,
      download_url = DOWNLOAD_URL,
      license = LICENSE,
      packages = find_packages(),
      include_package_data = True,
    #   classifiers = [
    #     'Development Status :: 4 - Beta',
    #     'Environment :: Console',
    #     'Intended Audience :: Science/Research',
    #     'Natural Language :: English',
    #     'Programming Language :: Python :: 3.7',
    #     'Programming Language :: Python :: 3 :: Only',
    #     'Topic :: Scientific/Engineering :: Astronomy'],
     )