from collections import defaultdict
debugoutput=defaultdict(lambda :False)
debugoutput['path']='dbout'
from .optimize import *
from .optutils import *
from .acquisitions import *
from .reccomenders import *
from .config import *
from .GPdc import *
from .ESutils import *
from .choosers import *
