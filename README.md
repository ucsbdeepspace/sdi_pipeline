#SDI Pipeline
UCSB’s astronomical data processing pipeline

## Installation
Currently only prepared for use on `starmaker`

## Usage
```
$ python3 -m /path/to/sdi /path/to/fits/ 

To save difference images, add "-s" and desired file path to a save folder: 

$ python3 -m /path/to/sdi -s /path/to/fits/ /path/to/save/folder/
```
Outputs a `srcs.pkl` file representing a catalog of transient candidates.
