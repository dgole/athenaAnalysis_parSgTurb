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
################################################################################
pathBase = str(sys.argv[1])
#keyList = ['drho', 'dvx', 'dvy', 'dvz']
keyList = ['drho']
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/acf/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
################################################################################
for key in keyList:
    acfMean = reader3d.acfMean(do3d, key)
    reader3d.plotAcf(acfMean, figNum=0)
    tools.saveAndClear(pathSave + key + '_stAvg_viridis.png',  figNum=0)
    reader3d.plotAcfDiverging(acfMean, figNum=0)
    tools.saveAndClear(pathSave + key + '_stAvg_coolwarm.png', figNum=0)
################################################################################




















#
