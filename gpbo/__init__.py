import os
from . import core
from . import test
from . import examples



with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'VERSION')) as version_file:
    __version__ = version_file.read().strip()

def about():
    return "I'm a function in __init__.py"