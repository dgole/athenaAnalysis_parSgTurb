#!/usr/bin/python
from __future__ import unicode_literals
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import os
import math
from scipy import fftpack as fft
import sys
import athenaReader as reader
import resource
from matplotlib.backends.backend_pdf import PdfPages
m.rcParams['text.usetex'] = True
m.rcParams['text.latex.unicode'] = True

# functions list
#__init__(self, path, iStart, iEnd, baseName="StratCooling", dt=0.1):
#addCol(self, funcName, label, *args, **kwargs):
#stPlot(self, col, cmapType="viridis", logOption=0, save=None, savePath=None, clim1=None, clim2=None):
#profile(self, col, nStart, nEnd, logOption=0, save=None, savePath=None, legendLabel=None):
#timeEvo(self, col, zStart, zEnd, logOption=0, save=None, savePath="", legendLabel=None):
#makeStandardPlotsPdf(self, t1, t2, t3):

# initialize 1d and 3d data object lists
run1d = []; run3d = [];


###############################################################
####################        1d        #########################
###############################################################

fileList = os.listdir("./out" + str(sys.argv[1])+ "/1d/")
nFiles = 0;
for i in range(len(fileList)):
	if (fileList[i][:12]=="StratCooling"): nFiles=nFiles+1;
run1d.append(reader.Data1d("./out" + str(sys.argv[1]) + "/1d/", 0, nFiles))
t1=nFiles/10/2; t2=nFiles/10;


run1d[0].makeStandardPlotsPdf(0, t1, t2)
run1d[0].stPlot(3, cmapType="viridis", save="pdf", clim1=0.4, clim2=10.0)
run1d[0].stPlot(17, cmapType="viridis", save="pdf")

# add column and plot it 
run1d[0].addCol(reader.dv1d, r"$ \delta v / c_s $")
run1d[0].profile(19, t1, t2, logOption=1, save="pdf")
run1d[0].stPlot(19, logOption=1, save="pdf")

run1d[0].addCol(reader.transEff, "Transport Effeciency")
run1d[0].profile(20, t1, t2, logOption=1, save="pdf")
run1d[0].stPlot(20, logOption=1, clim1=-2, save="pdf")

# multiple timeEvos plot
col=18
run1d[0].timeEvo(col, -3.0, -2.0, logOption=1, legendLabel="-z")
run1d[0].timeEvo(col, -0.5,  0.5, logOption=1, legendLabel="midplane")
run1d[0].timeEvo(col,  2.0,  3.0, logOption=1, legendLabel="+z")
plt.savefig(run1d[0].pdfName, format='pdf'); plt.clf();

# multiple profiles plot
run1d[0].profile(18, t1, t2, logOption=1, legendLabel="maxwell")
run1d[0].profile(10, t1, t2, logOption=1, legendLabel="reynolds")
plt.savefig(run1d[0].pdfName, format='pdf'); plt.clf();




###############################################################
####################        3d        #########################
###############################################################
'''
run3d.append(reader.Data3d("./out3020/3dComb/", 0, 100))
t1=70; t2=100;

run3d[0].addCol(reader.maxwell3d, r"$B_x B_y$")
run3d[0].profile(11, t1, t2, logOption=1, save="pdf")
run3d[0].stPlot(11, logOption=1, save="pdf", clim1=-7.0)
run3d[0].timeEvo(11, -3.0, -2.0, logOption=1, legendLabel="-z")
run3d[0].timeEvo(11, -0.5,  0.5, logOption=1, legendLabel="midplane")
run3d[0].timeEvo(11,  2.0,  3.0, logOption=1, legendLabel="+z")
plt.savefig(run3d[0].pdfName, format="pdf"); plt.clf();

#run3d[0].addCol(reader.dv3d, r"$\delta v$")
#run3d[0].profile(12, 70, 100, logOption=1, save="png"); plt.clf();
run3d[0].addCol(reader.zeros3d, "place holder")

#run3d[0].addCol(reader.drho3d, r"$\delta \rho / \rho$")
#run3d[0].profile(13, 0, 100, logOption=1, save="png"); plt.clf();
run3d[0].addCol(reader.zeros3d, "place holder")

run3d[0].addCol(reader.keFluxx, "KE Flux x")
run3d[0].addCol(reader.keFluxy, "KE Flux y")
run3d[0].addCol(reader.keFluxz, "KE Flux z")
run3d[0].profile(14, t1, t2, logOption=0, save="pdf")
run3d[0].profile(15, t1, t2, logOption=0, save="pdf")
run3d[0].profile(16, t1, t2, logOption=0, save="pdf")
run3d[0].stPlot(14, logOption=1, clim1=-7, save="pdf")
run3d[0].stPlot(15, logOption=1, clim1=-7, save="pdf")
run3d[0].stPlot(16, logOption=1, clim1=-7, save="pdf")

'''






# close all open pdfs
for obj in run1d: obj.pdfName.close()
#for obj in run3d: obj.pdfName.close()






################ 1d cols #################### 
#0 x3
#1 dens
#2 pressure
#3 temperature
#4 E
#5 Etot
#6 KEx
#7 KEy
#8 KEz
#9 KE
#10 Reynolds
#11 MEx
#12 MEy
#13 MEz
#14 ME
#15 Bx
#16 By
#17 Bz
#18 Maxwell

################ 3d cols #################### 
#0 x
#1 y
#2 z
#3 dens
#4 v1
#5 v2
#6 v3
#7 P
#8 B1
#9 B2
#10 B3













































