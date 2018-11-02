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
from matplotlib.backends.backend_pdf import PdfPages

doGas1d  = 1
doParAvg = 1

# data path and output path
path        = '../../data/' + str(sys.argv[1])
savePath    = path + 'plots/'
if not os.path.exists(savePath): os.makedirs(savePath)

def saveAndClear(plotName, figNum=0):
	plt.figure(figNum)
	plt.savefig(savePath + plotName, bbox_inches='tight')
	plt.clf()

# 1d GAS
if doGas1d == 1:
	dataPath1d  = path + '1d/'
	do1d = reader1d.Data1d(dataPath1d)
	do1d.addCol(reader1d.dv, 'dv', r'$\delta v$')
	for key in do1d.data.keys():
		reader1d.stPlot(do1d, key)
		saveAndClear("gas_ST_" + key + ".png")
		reader1d.profile(do1d, key)
		saveAndClear("gas_profile_" + key + ".png")
		reader1d.timeEvo(do1d, key)
		saveAndClear("gas_timeEvo_" + key + ".png")

# PARTICLE AVERAGES
if doParAvg == 1:
	doPhst = readerPhst.DataPhst(path)
	for key in ['zvar', 'vzvar']:
		readerPhst.timeEvo(doPhst, key)
		saveAndClear("par_timeEvo_" + key + ".png")






















#
