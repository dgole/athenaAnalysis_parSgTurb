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
pathBase    = str(sys.argv[1])
runNameList = ['run220/bin', 'run221/bin', 'run222/bin', 'run223/bin', 'run224/bin']
labelList   = ['64_1_32', '128_1_32', '64_1_64', '128_1_64', '128_1_128']
colorList   = ['b', 'r', 'b', 'r', 'b']
sytleList   = ['-', '-', '--', '--', ':']
pathSave = pathBase + 'plots/multiPlots3d_220/'
################################################################################
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3dList    = []
for n in range(len(runNameList)):
	path3d     = pathBase + runNameList[n] + '/3d/'
	do3dList.append(reader3d.Data3d(path3d))
################################################################################
key = 'dv'
for n in range(len(runNameList)):
    reader3d.profile(do3dList[n], key, figNum=0, absAvg=1, absPlot=1,
    legendLabel = labelList[n], color=colorList[n], linestyle=styleList[n])
plt.legend()
plt.ylim(1.e-3,1.e0)
tools.saveAndClear(pathSave + 'profile_' + key + '.png', figNum=0)



#
