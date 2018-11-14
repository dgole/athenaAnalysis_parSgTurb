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
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/quickAcf/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
#############################################################
shiftMax   = do3d.nz // 2
sCutFactor = do3d.nz // 8
tCutFactor = 6
ziList = range(sCutFactor//2, do3d.nz, sCutFactor)
nList  = range(do3d.nt-6*4-1, do3d.nt, tCutFactor)
for key in keyList:
    acfMean = reader3d.acfMean(do3d, key, shiftMax, ziList, nList)
    reader3d.plotAcf(acfMean, figNum=0)
    tools.saveAndClear(pathSave + key + '_stAvg.png', figNum=0)
#############################################################




















#
