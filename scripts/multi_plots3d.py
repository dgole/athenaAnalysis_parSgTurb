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
#pathBase    = str(sys.argv[1])
#runNameList = ['run120',  'run121', 'run122', 'run123']
#labelList   = ['64_1_32', '128_1_32', '128_1_64', '256_1_128']
#colorList   = ['r','g','b','m']
#styleList   = ['-', '-', '-', '-']
#pathSave    = pathBase + 'plots/multiPlots3d_120/'
################################################################################
#pathBase    = '../../data/kspaceTest/'
#runNameList = ['run150', 'run151', 'run152', 'run153', 'run154']
#labelList   = ['1/3', '-2/3', '-5/3', '-8/3', '-11/3']
#colorList   = ['r', 'g', 'm', 'b', 'k']
#pathSave = pathBase + 'plots/multiPlots3d_150/'
#styleList   = ['-', '-', '-', '-', '-']
################################################################################
pathBase    = '../../data/kspaceTest/'
runNameList = ['run160', 'run161', 'run162', 'run163', 'run164']
labelList   = ['dedt=1.e-8', 'dedt=1.e-7', 'dedt=1.e-6', 'dedt=1.e-5', 'dedt=1.e-4']
colorList   = ['r', 'g', 'm', 'b', 'k']
pathSave = pathBase + 'plots/multiPlots3d_160/'
styleList   = ['-', '-', '-', '-', '-']
################################################################################
nStart      = 1
nEnd        = 5
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3dList    = []
for n in range(len(runNameList)):
	path3d     = pathBase + runNameList[n] + '/3d/'
	do3dList.append(reader3d.Data3d(path3d))
################################################################################
key = 'dv'
for n in range(len(runNameList)):
    reader3d.profile(do3dList[n], key, figNum=0, absAvg=1, absPlot=1,
    legendLabel = labelList[n], color=colorList[n], linestyle=styleList[n],
	nStart=nStart, nEnd=nEnd)
plt.legend()
plt.ylim(1.e-3,1.e0)
tools.saveAndClear(pathSave + 'profile_' + key + '.png', figNum=0)



#
