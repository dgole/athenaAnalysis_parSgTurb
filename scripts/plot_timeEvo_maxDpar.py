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
pathBase = str(sys.argv[1])
nStart   = int(sys.argv[2])
nEnd     = int(sys.argv[3])
alpha    = float(sys.argv[4])
dvBase   = np.sqrt(alpha)
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/plots3d/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
################################################################################
# max dpar
plt.figure(0)
print("making max dpar plot...")
sys.stdout.flush()
key = 'dpar'
plotDataList = []
for n in range(do3d.nt):
	print("reading in time step " + str(n) + "...")
	sys.stdout.flush()
	data = do3d.get3d(key, n)
	plotDataList.append(np.amax(data))
plotData = np.asarray(plotDataList)
print(plotData)
plt.xlabel(r'$t \Omega$')
plt.ylabel('MAX(' + do3d.header[key] + ')')
plt.semilogy(do3d.t, plotData, color='k')
plt.axvline(20, color='gray', lineStyle=':')
plt.axhline(plotData[20], color='gray', lineStyle=':')
plt.xlim(0.0, np.amax(do3d.t))
tools.saveAndClear(pathSave + 'timeEvo_maxDPar.png', figNum=0)





#
