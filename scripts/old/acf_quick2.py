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
#############################################################
pathBase = str(sys.argv[1])
if len(sys.argv) == 2:
    keyList = [
              'drho'#, 'dvx', 'dvy', 'dvz', 'dv'
              #,'drhoNorm', 'dvxNorm', 'dvyNorm', 'dvzNorm', 'dvNorm'
              ]
else:
    keyList = []
    for arg in sys.argv[5:]: keyList.append(str(arg))
#############################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/quickAcf2/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
#############################################################
key = 'drho'; n=100; extent=[-0.1,0.1,-0.1,0.1];
nn = do3d.nx//2
acf3d = reader3d.acf3d(do3d, key, n)
acf2d = np.mean(acf3d, axis=2)
plotData  = np.transpose(np.fliplr(acf2d))
plotData2 = np.zeros_like(plotData)
plotData2[:nn,:nn] = plotData[nn:,nn:] # 3 into 1
plotData2[:nn,nn:] = plotData[nn:,:nn] # 4 into 2
plotData2[nn:,nn:] = plotData[:nn,:nn] # 1 into 3
plotData2[nn:,:nn] = plotData[:nn,nn:] # 2 into 4
plt.figure(0)
plt.imshow(plotData2, extent=extent, aspect=1.0, cmap=plt.get_cmap('coolwarm'))
plt.colorbar()
#plt.clim(-1,1)
plt.tight_layout()
tools.saveAndClear(pathSave + key + '_stAvg.png', figNum=0)











#############################################################




















#
