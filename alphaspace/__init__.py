
try:
    import numpy as np
except:
    print("Numpy is needed")
    exit()

try:
    import scipy
except:
    print("scipy is needed")
    exit()
try:
    import mdtraj
except:
    print("mdtraj is needed")
    exit()

try:
    import nglview as nv
except:
    print("nglview not found, jupyter notebook visualization will be disabled")

from AS_Universe import AS_Universe
from AS_Funct import *

