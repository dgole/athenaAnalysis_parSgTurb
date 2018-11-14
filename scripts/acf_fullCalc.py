#!/usr/bin/python
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import math
import sys
sys.path.append('../python')
import athenaReader3d as reader3d
import athenaTools as tools
from matplotlib.backends.backend_pdf import PdfPages
#############################################################
pathBase = str(sys.argv[1])
if len(sys.argv) == 2:
    keyList = [
              'drho', 'dvx', 'dvy', 'dvz', 'dv',
              'drhoNorm', 'dvxNorm', 'dvyNorm', 'dvzNorm', 'dvNorm'
              ]
else:
    keyList = []
    for arg in sys.argv[2:]: keyList.append(str(arg))
#############################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'acf4d/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
#############################################################
shiftMax = do3d.nx // 2
for key in keyList:
    acf4d = reader3d.acf4d(do3d, key, shiftMax)
    np.save(pathSave + key + '.npy', acf4d)
#############################################################




















#
