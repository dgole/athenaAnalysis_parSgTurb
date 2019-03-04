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
pathBase = '../../data/turbTest/'
runNameList = ['run36', 'run35', 'run30', 'run31', 'run32', 'run33', 'run34']
alphaInList = [1.e0,     1.e-1,   1.e-2,   1.e-3,   1.e-4,   1.e-5,   1.e-6 ]
colorList   = ['b',  'tab:orange', 'g',    'r',     'tab:purple', 'tab:brown', 'tab:pink'  ]
pathSave = pathBase + 'plots/3danalysis30/'

pathBase = '../../data/turbTest/'
runNameList = ['run30', 'run34']
alphaInList = [1.e-2,   1.e-6 ]
tsList      = [1.e-1,   1.e-1 ]
colorList   = ['b',  'tab:orange', 'g',    'r',     'tab:purple', 'tab:brown', 'tab:pink'  ]
pathSave = pathBase + 'plots/3danalysis30/'
################################################################################
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3dList    = []
for n in range(len(runNameList)):
	path3d     = pathBase + runNameList[n] + '/3d/'
	do3dList.append(reader3d.Data3d(path3d))
################################################################################
# simple averages of quantities
key = 'vz'
for n in range(len(runNameList)):
    reader3d.profile(do3dList[n], key, figNum=0, absAvg=0, absPlot=0,
    legendLabel = r'$\alpha_{in}=$'+str(alphaInList[n]), color=colorList[n])
plt.legend()
plt.ylim(-1.e-2,1.e-2)
tools.saveAndClear(pathSave + 'profileRealAvg_' + key + '.png', figNum=0)
