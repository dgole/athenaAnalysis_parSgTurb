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
import athenaReaderIsothermal as reader
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




###############################################################
####################        1d        #########################
###############################################################
run1d = [];
run1d.append(reader.Data1d("/home/dan/ResearchLocal/prod_runs/athenaJanus702emfFix/1d/", 0, 5000))
run1d[0].stPlot(1, cmapType="viridis", logOption=1, save="pdf")

run1d[0].addCol(reader.vx, r"$V_x$")
run1d[0].addCol(reader.vy, r"$V_y$")
run1d[0].addCol(reader.vz, r"$V_z$")

run1d[0].stPlot(16, cmapType="viridis", logOption=1, save="pdf", clim1=-3.0)









for obj in run1d: obj.pdfName.close()




################ 1d cols #################### 
#0 x3
#1 dens
#2 pressure
#3 KEx
#4 KEy
#5 KEz
#6 KE
#7 Reynolds
#8 MEx
#9 MEy
#10 MEz
#11 ME
#12 Bx
#13 By
#14 Bz
#15 Maxwell
#16 Vx
#17 Vy
#18 Vz






































