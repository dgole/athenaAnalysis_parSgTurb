#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
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
pathAcf  = pathBase + 'acf4d/'
pathSave = pathBase + 'plots/fullAcf/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
#############################################################
for key in keyList:
    fileName = pathAcf + key + '.npy'
    acf4d   = np.load(fileName)
    acfMean = np.mean(acf4d[do.nt//2:], axis=(0,1))
    reader3d.plotAcf(acfMean, figNum=0)
    tools.saveAndClear(pathSave + key + '_stAvg.png', figNum=0)
#############################################################




















#
