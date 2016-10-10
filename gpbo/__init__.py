import os
from . import core
from . import test
from . import examples
from . import _configs
from core import search


with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'VERSION')) as version_file:
    __version__ = version_file.read().strip()

import logging
logging.basicConfig(level=logging.DEBUG)

def about():
    return "I'm a function in __init__.py"

