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
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/plots3d/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
################################################################################
# simple averages of quantities
for key in ['vx', 'vy', 'vz', 'v']:
    reader3d.profile(do3d, key, figNum=0, absAvg=0, absPlot=0)
    tools.saveAndClear(pathSave + 'profileRealAvg_' + key + '.png', figNum=0)
################################################################################
# abs averages of quantities
for key in ['vx', 'vy', 'vz', 'v']:
    reader3d.profile(do3d, key, figNum=0, absAvg=1, absPlot=1)
    tools.saveAndClear(pathSave + 'profileAbsAvg_' + key + '.png', figNum=0)
################################################################################
# perts
for key in ['drho', 'dvx', 'dvy', 'dvz', 'dv']:
    reader3d.profile(do3d, key, figNum=0)
    tools.saveAndClear(pathSave + 'profilePert_' + key + '.png', figNum=0)
################################################################################
# normalized perts
for key in ['drhoNorm', 'dvxNorm', 'dvyNorm', 'dvzNorm', 'dvNorm']:
    reader3d.profile(do3d, key, figNum=0)
    tools.saveAndClear(pathSave + 'profilePertNorm_' + key + '.png', figNum=0)
################################################################################




'''
for key in ['dvx','dvy','dvz', 'dv']:
    reader3d.profile(do3d, key, figNum=0, legendLabel=do3d.header[key])
plt.legend()
plt.ylim(1.e-3, 1.e-1)
tools.saveAndClear(pathSave + 'multiVelocity_.png', figNum=0)
################################################################################
for key in ['vx', 'vy', 'vz']:
    reader3d.profile(do3d, key, figNum=0, absAvg=1, absPlot=1)
    tools.saveAndClear(pathSave + 'profileAbs_' + key + '.png', figNum=0)
################################################################################
for key in ['vx', 'vy', 'vz']:
    reader3d.profile(do3d, key, figNum=0, absAvg=0, absPlot=0)
    tools.saveAndClear(pathSave + 'profile_' + key + '.png', figNum=0)
################################################################################
'''



















#
