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
kStart   = int(sys.argv[4])
kEnd     = int(sys.argv[5])
pathIn   = pathBase + 'pspecData/'
pathSave = pathBase + 'plots/pspec/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
################################################################################
vExpo = -1.833
eExpo = 2.0*(vExpo)+2.0
print("velocity spectrum PL exponent is " + str(vExpo))
print("KE power spectrum PL exponent is " + str(eExpo))
################################################################################
fileName = "pspecData_" + str(nStart)+"_"+str(nEnd)+"_"+str(kStart)+"_"+str(kEnd)+".npy"
inData = np.load(pathIn + fileName)
freqs  = inData[0,:]
psk  = inData[1,:]
psk *= np.power(freqs, -eExpo)
psk /= np.mean(psk)
print(freqs)
print(psk)
plt.loglog(freqs, psk)
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.ylim(1.e-2,1.e2)
plt.xlim(freqs[1],freqs[-1])
tools.saveAndClear(pathSave + 'adjustedPspec_2D_' +
                   str(kStart) + '_' + str(kEnd)  +
                   '.png', figNum=0)









#
