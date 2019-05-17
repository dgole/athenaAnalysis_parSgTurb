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
pathBase    = '../../data/kspaceTest/'
pathSave = pathBase + 'plots/'
dvDesired = 3.2e-2
plt.figure(0)
################################################################################
#dedt = np.asarray([1.e-8,1.e-7,1.e-6,1.e-5,1.e-4])
#dv   = np.asarray([0.001244,0.003755,0.01085,0.02836,0.06710])

#dedt    = np.asarray([1.e-7,   1.e-6,  1.e-5,  1.e-4])
#dv      = np.asarray([0.00659, 0.0168, 0.0391, 0.0864])

#dedt    = np.asarray([1.e-7,  1.e-6,  1.e-5,  1.e-4, 1.e-3])
#dv      = np.asarray([0.00126, 0.00374, 0.0107, 0.0285, 0.071])

dedt    = np.asarray([1.e-7,  1.e-6,  1.e-5,  1.e-4])
dv      = np.asarray([0.00185, 0.00549, 0.0154, 0.0399])

dedtLog = np.log10(dedt)
dvLog   = np.log10(dv)

plt.loglog(dedt, dv, 'ko', markersize=5)
plt.loglog(dedt, dv, color=(0,0,0,0.2))

dvDesiredLog   = np.log10(dvDesired)
dedtDesiredLog = np.interp(dvDesiredLog, dvLog, dedtLog)
dedtDesired    = np.power(10.0,dedtDesiredLog)
plt.loglog(dedtDesired, dvDesired, 'ro', markersize=5)
plt.axvline(x=dedtDesired, linestyle='--', color='r')
plt.axhline(y=dvDesired, linestyle='--', color='r')

plt.title('drive at dedt=' + str(dedtDesired) +
		 ' to achieve dv=' + str(dvDesired))
plt.xlabel('dedt')
plt.ylabel('total dv')
tools.saveAndClear(pathSave + 'alphaPred.png', figNum=0)
################################################################################
for alphaDesired in [1.e-5, 3.2e-5, 1.e-4, 3.2e-4, 1.e-3, 3.2e-3]:
	dvDesired      = np.sqrt(alphaDesired)
	dvDesiredLog   = np.log10(dvDesired)
	dedtDesiredLog = np.interp(dvDesiredLog, dvLog, dedtLog)
	dedtDesired    = np.power(10.0,dedtDesiredLog)
	print('alpha='+str(alphaDesired) + ' --> ' +
	      'dv=' + str(np.round(dvDesired,5)) + ' --> ' +
		  'dedt=' + str(dedtDesired))










#
