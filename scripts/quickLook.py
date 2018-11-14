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
pathSave = pathBase + 'plots/quickLook/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
#########################################################
# 1d GAS
if doGas1d == 1:
	do1d = reader1d.Data1d(path1d)
	do1d.addCol(reader1d.dv, 'dv', r'$\delta v$')
	do1d.addCol(reader1d.dvz, 'dvz', r'$\delta v_z$')
	for key in do1d.data.keys():
		reader1d.stPlot(do1d, key)
		tools.saveAndClear(pathSave + "gas_ST_" + key + ".png")
		reader1d.profile(do1d, key)
		tools.saveAndClear(pathSave + "gas_profile_" + key + ".png")
		reader1d.timeEvo(do1d, key)
		tools.saveAndClear(pathSave + "gas_timeEvo_" + key + ".png")
# PARTICLE AVERAGES
if doParAvg == 1:
	doPhst = readerPhst.DataPhst(pathBase)
	for key in ['zvar', 'vzvar']:
		readerPhst.timeEvo(doPhst, key)
		tools.saveAndClear(pathSave + "par_timeEvo_" + key + ".png")






















#
