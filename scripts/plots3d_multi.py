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
#pathBase    = '../../data/prodRuns/'
#runNameList = ['run301_noPar', 'run302_noPar', 'run303_noPar']
#abelList   = [r'$\alpha_{in}=10^{-3}$',
	           #r'$\alpha_{in}=10^{-3.5}$',
			   #r'$\alpha_{in}=10^{-4}$']
#colorList   = ['r', 'g', 'm', 'b']
#pathSave = pathBase + 'plots/multiPlots3d_300_noPar/'
#alphaDesiredList = [1.e-3, 3.2e-4, 1.e-4]
################################################################################
pathBase    = '../../data/kspaceTest/'
runNameList = ['run410', 'run411', 'run412', 'run413']
labelList   = [r'$\alpha_{in}=10^{-3}$',
	           r'$\alpha_{in}=10^{-3.5}$',
	           r'$\alpha_{in}=10^{-3.5}$',
			   r'$\alpha_{in}=10^{-4}$']
colorList   = ['r', 'g', 'm', 'b']
pathSave = pathBase + 'plots/multiPlots3d_400/'
alphaDesiredList = [1.e-3, 3.2e-4, 1.e-4, 1.e-4]
################################################################################
nStart      = 2
nEnd        = 7
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3dList    = []
for n in range(len(runNameList)):
	path3d     = pathBase + runNameList[n] + '/3d/'
	do3dList.append(reader3d.Data3d(path3d))
################################################################################
key = 'dv'
for n in range(len(runNameList)):
	dvMean = reader3d.profile(do3dList[n], key, figNum=0, absAvg=1, absPlot=1,
	legendLabel = labelList[n],	color=colorList[n], linestyle='-',
	nStart=nStart, nEnd=nEnd)
	plt.axhline(y=np.sqrt(alphaDesiredList[n]), lineStyle='--', color=colorList[n])
	alphaObs = np.square(dvMean)
	print(r'$\alpha_{obs}=$' + str(alphaObs))
plt.legend()
plt.ylim(3.e-3,6.e-2)
tools.saveAndClear(pathSave + 'profile_' + key + '.png', figNum=0)



#
