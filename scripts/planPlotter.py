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
import planOutputReader as readerPlan
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
################################################################################
scatters = False
################################################################################
# paths
pathBase = str(sys.argv[1])
pathPlan = pathBase + 'planOutput/'
pathSave = pathBase + 'plots/plan/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
################################################################################
# 1    2    3      4    5    6
# npar mass r_hill xcom ycom zcom
doPlan = readerPlan.DataPlan(pathPlan)
################################################################################
# n clumps over time
plt.figure(0)
readerPlan.nClumpsTimeEvo(doPlan, figNum=0)
tools.saveAndClear(pathSave + "nClumps_timeEvo.png", figNum=0)
################################################################################
# mass cumulative hist
readerPlan.plotCumMassHist(doPlan)
tools.saveAndClear(pathSave + "hist_cumulative.png", figNum=0)
################################################################################
# mass differential hist
readerPlan.plotDiffMassHist(doPlan)
tools.saveAndClear(pathSave + "hist_differential.png", figNum=0)
################################################################################
# fraction of particle mass in planetesimals
readerPlan.planMassFracTimeEvo(doPlan)
tools.saveAndClear(pathSave + "massfrac_timeEvo.png", figNum=0)
################################################################################
# scatter plot of particle locations
if scatters:
	plt.figure(0)
	for n in range(doPlan.nFirstClump, doPlan.nt):
		readerPlan.scatterPlotXZ(doPlan, n)
		tools.saveAndClear(pathSave + 'scatter_2d_XZ' + '_' + str(n) + '.png', figNum=0, dpi=100)
		readerPlan.scatterPlotXY(doPlan, n)
		tools.saveAndClear(pathSave + 'scatter_2d_XY' + '_' + str(n) + '.png', figNum=0, dpi=100)
		readerPlan.scatterPlotXYZ(doPlan, n)
		tools.saveAndClear(pathSave + 'scatter_3d' + '_' + str(n) + '.png', figNum=0, dpi=100, bboxOption=0)
################################################################################
# calculate p value (power law index)
plt.figure(0)
readerPlan.pValuePlot(doPlan)
tools.saveAndClear(pathSave + "powerLawIndex.png", figNum=0)
################################################################################
















#
