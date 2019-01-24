#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import athenaReader1d as reader1d
import athenaReaderPhst as readerPhst
import athenaTools as tools
from matplotlib.backends.backend_pdf import PdfPages
#########################################################
doGas1d  = 1
doParAvg = 1
# paths
pathBase = str(sys.argv[1])
path1d   = pathBase + '1d/'
pathSave = pathBase + 'plots/sloshChecker/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
#########################################################
# 1d GAS
if doGas1d == 1:
	do1d = reader1d.Data1d(path1d)
	for tStart in np.arange(0, 2, 0.1):
		reader1d.profile(do1d, 'rho', tStart=tStart, tEnd=tStart+0.1, legendLabel=tStart)
	plt.legend()
	tools.saveAndClear(pathSave + "profiles_rho.png")
	reader1d.timeEvo(do1d, 'rho', z1=0.0, z2=0.01)
	tools.saveAndClear(pathSave + "timeEvo_rho.png")





















#
