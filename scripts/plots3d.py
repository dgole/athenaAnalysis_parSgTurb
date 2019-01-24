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
for key in ['vx', 'vy', 'vz']:
    reader3d.profile(do3d, key, figNum=0, absAvg=0, absPlot=0)
    tools.saveAndClear(pathSave + 'profileRealAvg_' + key + '.png', figNum=0)
    reader3d.timeEvo(do3d, key, figNum=0, absAvg=0, absPlot=0)
    tools.saveAndClear(pathSave + 'timvoRealAvg_' + key + '.png', figNum=0)
################################################################################
# abs averages of quantities
for key in ['vx', 'vy', 'vz']:
    reader3d.profile(do3d, key, figNum=0)
    tools.saveAndClear(pathSave + 'profileAbsAvg_' + key + '.png', figNum=0)
    reader3d.timeEvo(do3d, key, figNum=0)
    tools.saveAndClear(pathSave + 'timeEvoAbsAvg_' + key + '.png', figNum=0)
################################################################################
# perts
for key in ['drho', 'dvx', 'dvy', 'dvz']:
    reader3d.profile(do3d, key, figNum=0)
    tools.saveAndClear(pathSave + 'profilePert_' + key + '.png', figNum=0)
    reader3d.timeEvo(do3d, key, figNum=0)
    tools.saveAndClear(pathSave + 'timeEvoPert_' + key + '.png', figNum=0)
################################################################################
# normalized perts
for key in ['drhoNorm', 'dvxNorm', 'dvyNorm', 'dvzNorm', 'dvNorm']:
    reader3d.profile(do3d, key, figNum=0)
    tools.saveAndClear(pathSave + 'profilePertNorm_' + key + '.png', figNum=0)
    reader3d.timeEvo(do3d, key, figNum=0)
    tools.saveAndClear(pathSave + 'timeEvoPertNorm_' + key + '.png', figNum=0)
################################################################################

















#
