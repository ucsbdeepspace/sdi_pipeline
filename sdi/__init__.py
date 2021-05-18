from .align import align
from .combine import combine
from .display import display
from .fitsio import read, write
from .subtract import subtract
from .extract import extract
from .ref import ref
from . import test

# in order to test click commands
from .fitsio import read_cmd as _read_cmd
from .fitsio import write_cmd as _write_cmd
from .align import align_cmd as _align_cmd
from .combine import combine_cmd as combine_cmd
from .extract import extract_cmd as _extract_cmd
from .subtract import subtract_cmd as _subtract_cmd
from . import test_cmd

import os
_scidir = os.path.join(os.path.dirname(__file__), "test/fixtures/science")
_resdir = os.path.join(os.path.dirname(__file__), "test/fixtures/residuals")
