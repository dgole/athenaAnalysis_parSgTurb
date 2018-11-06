#!/usr/bin/python
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import athenaReader1d as reader1d
import athenaReaderPhst as readerPhst
import athenaTools as tools
from matplotlib.backends.backend_pdf import PdfPages
####################################################

pathBase = '../../data/turbTest/'
pathSave = pathBase + 'plots/turbAnalysis/'
if not os.path.exists(pathSave): os.makedirs(pathSave)

runNameList = ['run10', 'run11', 'run12', 'run13', 'run14']
alphaInList = [1.e-2,    1.e-3,   1.e-4,   1.e-5,   1.e-6]
colorList   = ['b',  'tab:orange', 'g',    'r',     'tab:purple']

# read 1d and phst files in for all runs
do1dList    = []
doPhstList  = []
for n in range(len(runNameList)):
	path1d     = pathBase + runNameList[n] + '/1d/'
	pathParHst = pathBase + runNameList[n] + '/'
	do1dList  .append(reader1d.Data1d(path1d))
	doPhstList.append(readerPhst.DataPhst(pathParHst))
	do1dList[n].addCol(reader1d.dv,     'dv',     r'$\delta v$')
	do1dList[n].addCol(reader1d.dvz,    'dvz',    r'$\delta v_z$')
	do1dList[n].addCol(reader1d.alphaz, 'alphaz', r'$\delta v_z$')

# plot dv vs expected value base on alpha_in
key = 'dv'
for n in range(len(do1dList)):
	do      = do1dList[n]
	alphaIn = alphaInList[n]
	color   = colorList[n]
	reader1d.timeEvo(do, key, legendLabel=r'$\alpha_{in}=$'+str(alphaIn), logForce=1)
	plt.axhline(y=np.sqrt(alphaIn), linestyle='--', color=color)
plt.legend()
tools.saveAndClear(pathSave + "gas_timeEvo_" + key + ".png")

# plot particle dz vs expected value based on alpha ~ dvz^2
for n in range(len(do1dList)):
	do         = do1dList[n]
	alphaIn    = alphaInList[n]
	color      = colorList[n]
	doPhst     = doPhstList[n]
	alphazMean = np.mean(do.data['alphaz'][do.nt//2:])
	parh       = np.sqrt(alphazMean/0.1)
	readerPhst.timeEvo(doPhst, 'zvar', legendLabel=r'$\alpha_{in}=$'+str(alphaIn), logForce=1)
	plt.axhline(y=parh, linestyle='--', color=color)
plt.axhline(y=do.dz, linestyle='--', color='k')
plt.legend()
tools.saveAndClear(pathSave + "par_" + 'scaleHeightComparison' + ".png")

# plot particle dz vs expected value based on alpha ~ dvz^2
for n in range(len(do1dList)):
	do         = do1dList[n]
	alphaIn    = alphaInList[n]
	color      = colorList[n]
	doPhst     = doPhstList[n]
	alphazMean = np.mean(do.data['alphaz'][do.nt//2:, do.nz//2-1:do.nz//2+2])
	parh       = np.sqrt(alphazMean/0.1)
	readerPhst.timeEvo(doPhst, 'zvar', legendLabel=r'$\alpha_{in}=$'+str(alphaIn), logForce=1)
	plt.axhline(y=parh, linestyle='--', color=color)
plt.axhline(y=do.dz, linestyle='--', color='k')
plt.legend()
tools.saveAndClear(pathSave + "par_" + 'scaleHeightComparison2' + ".png")




# PARTICLE AVERAGES
#doPhst = readerPhst.DataPhst(pathBase)
#for key in ['zvar', 'vzvar']:
	#readerPhst.timeEvo(doPhst, key)
	#tools.saveAndClear(pathSave + "par_timeEvo_" + key + ".png")






















#
