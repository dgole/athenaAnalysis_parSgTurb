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
if str(sys.argv[1])=='11?':
	pathBase    = str(sys.argv[2]) 
	runNameList = ['run110/bin','run111/bin','run112/bin']
	labelList   = ['64','128','256']
	tsList      = [0.3, 0.3, 0.3]
	colorList   = ['tab:blue', 'tab:orange', 'g', 'r']
	pathSave = pathBase + 'plots/turbAnalysis110/'
################################################################################
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3dList    = []
for n in range(len(runNameList)):
	path3d     = pathBase + runNameList[n] + '/3d/'
	do3dList.append(reader3d.Data3d(path3d))
################################################################################
key = 'dvx'
for n in range(len(runNameList)):
    reader3d.profile(do3dList[n], key, figNum=0, absAvg=1, absPlot=1,
    legendLabel = labelList[n], color=colorList[n])
plt.legend()
plt.ylim(1.e-3,1.e0)
tools.saveAndClear(pathSave + 'profile_' + key + '.png', figNum=0)
################################################################################
key = 'dvy'
for n in range(len(runNameList)):
    reader3d.profile(do3dList[n], key, figNum=0, absAvg=1, absPlot=1,
    legendLabel = labelList[n], color=colorList[n])
plt.legend()
plt.ylim(1.e-3,1.e0)
tools.saveAndClear(pathSave + 'profile_' + key + '.png', figNum=0)
################################################################################
key = 'dvz'
for n in range(len(runNameList)):
    reader3d.profile(do3dList[n], key, figNum=0, absAvg=1, absPlot=1,
    legendLabel = labelList[n], color=colorList[n])
plt.legend()
plt.ylim(1.e-3,1.e-1)
tools.saveAndClear(pathSave + 'profile_' + key + '.png', figNum=0)
################################################################################
key = 'dv'
for n in range(len(runNameList)):
    reader3d.profile(do3dList[n], key, figNum=0, absAvg=1, absPlot=1,
    legendLabel = labelList[n], color=colorList[n])
plt.legend()
plt.ylim(1.e-3,1.e0)
tools.saveAndClear(pathSave + 'profile_' + key + '.png', figNum=0)
################################################################################
key = 'dv'
for n in range(len(runNameList)):
    reader3d.timeEvo(do3dList[n], key, figNum=0, absAvg=1, absPlot=1,
    legendLabel = labelList[n], color=colorList[n])
plt.legend()
plt.ylim(1.e-3,1.e0)
tools.saveAndClear(pathSave + 'timeEvo_' + key + '.png', figNum=0)
