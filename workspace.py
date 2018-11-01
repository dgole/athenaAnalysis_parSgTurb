#!/usr/bin/python
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import os
import math
import sys
import athenaReader1d as reader1d
from matplotlib.backends.backend_pdf import PdfPages


path     = '../data/turbTest/run1/'
dataPath = path + '1d/'
savePath  = path + 'plots/'
if not os.path.exists(savePath): os.makedirs(savePath)
do1d    = reader1d.Data(dataPath)

for col in range(1,7):
	reader1d.stPlot(do1d, col)
	plt.savefig(savePath + "ST_" + str(col) + ".png", bbox_inches='tight')
	plt.clf()

col = 7
reader1d.stPlot(do1d, col, cmapType='coolwarm')
plt.savefig(savePath + "ST_" + str(col) + ".png", bbox_inches='tight')
plt.clf()

































# add column and plot it
#run1d[0].addCol(reader.dv1d, r"$ \delta v / c_s $")
#run1d[0].profile(19, t1, t2, logOption=1, save="pdf")
#run1d[0].stPlot(19, logOption=1, save="pdf")

#run1d[0].addCol(reader.transEff, "Transport Effeciency")
#run1d[0].profile(20, t1, t2, logOption=1, save="pdf")
#run1d[0].stPlot(20, logOption=1, clim1=-2, save="pdf")

# multiple timeEvos plot
#col=18
#run1d[0].timeEvo(col, -3.0, -2.0, logOption=1, legendLabel="-z")
#run1d[0].timeEvo(col, -0.5,  0.5, logOption=1, legendLabel="midplane")
#run1d[0].timeEvo(col,  2.0,  3.0, logOption=1, legendLabel="+z")
#plt.savefig(run1d[0].pdfName, format='pdf'); plt.clf();

# multiple profiles plot
#run1d[0].profile(18, t1, t2, logOption=1, legendLabel="maxwell")
#run1d[0].profile(10, t1, t2, logOption=1, legendLabel="reynolds")
#plt.savefig(run1d[0].pdfName, format='pdf'); plt.clf();
