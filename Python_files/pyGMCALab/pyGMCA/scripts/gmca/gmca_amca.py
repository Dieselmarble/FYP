# Useful packages

import sys
PYGMCALAB_PATH = "/Users/kevin/Documents/FYP/Python_files/pyGMCALab/pyGMCA"
sys.path.insert(1,PYGMCALAB_PATH)
from common import utils as bu
from bss.amca import amca
from bss.gmca import gmca
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}
plt.rcParams["figure.figsize"] = (20,10)
plt.rc('font', **font)

N_MC = 10 # Number of Monte-Carlo simulations