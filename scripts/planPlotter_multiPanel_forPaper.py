#!/usr/bin/python
import numpy as np
import matplotlib as m
#m.use('Agg')
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
import numpy.polynomial.polynomial as poly
################################################################################
# paths and CL args
pathBase = str(sys.argv[1])
nStart   = int(sys.argv[2])
nTot     = int(sys.argv[3])
xLim1    = (nStart-200)/10.0
t_vline1 = 16.7
t_vline2 = 18.8
pathPlan = pathBase + 'planOutput2/'
pathSave = pathBase + 'plots/planForPaper/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
################################################################################
# set up data object and other params
# 1    2    3      4    5    6
# npar mass r_hill xcom ycom zcom
doPlan1 = readerPlan.DataPlan(pathPlan, nStart=nStart, nTot=nTot, nPar=512*512*512, dt=0.1)
p_mle_master   = []
err_mle_master = []
################################################################################
def calcStuff(doPlan, n):
	# calculate the PL slope with the MLE and fit
	#doPlot = True
	try:
		mp, dndmp = readerPlan.getDiffMassHist(doPlan, n)
		p_mle, err_mle        = readerPlan.get_p_mle(doPlan, n)
		p_mle_master.append(p_mle)
		err_mle_master.append(err_mle)
	except:
		#doPlot = False
		p_mle_master.append(0.0)
		err_mle_master.append(0.0)

for n in range(nStart, doPlan1.nTot):
	calcStuff(doPlan1, n)

print(p_mle_master)
# actually make the plot
# Setup
fac = 0.45
fig = plt.figure(figsize=(19*fac, 15*fac), dpi=120)
ax = []
ax.append(plt.subplot2grid((4, 1), (0, 0), colspan=1, rowspan=2)) # 0 3 bottom plots
ax.append(plt.subplot2grid((4, 1), (2, 0), colspan=1)) # 1
ax.append(plt.subplot2grid((4, 1), (3, 0), colspan=1)) # 2

# plot p values over time
# p_mle
axNum = 0
ax[axNum].plot(doPlan1.time[nStart:], p_mle_master, 'k', linewidth=1)
pArr = np.asarray(p_mle_master); errArr = np.asarray(err_mle_master)
ax[axNum].fill_between(doPlan1.time[nStart:], pArr-errArr, pArr+errArr, color='gray', alpha=0.3)
# other stuff
ax[axNum].axhline(1.4, color=(0,0,0,0.2), linestyle='--')
ax[axNum].axhline(1.7, color=(0,0,0,0.2), linestyle='--')
ax[axNum].set_ylim(1.0,3.0)
ax[axNum].set_xlim(xLim1, doPlan1.tMax)
ax[axNum].set_ylabel("SPL Slope")
ax[axNum].axvline(t_vline1, color=(0,0,1,1), linestyle='--')
ax[axNum].axvline(t_vline2, color=(1,0,0,1), linestyle='--')

# plot mass frac in planetesimals over time
axNum = 1
mFracList = []
for item in doPlan1.peakArrayList[nStart:]:
	try:
		mNow = np.sum(item[:,2])
	except:
		mNow = 0.0
	mFracNow = mNow / doPlan1.mParTot
	mFracList.append(mFracNow)
ax[axNum].plot(doPlan1.time[nStart:], mFracList, 'k', linewidth=2)
ax[axNum].set_ylabel(r'$M_{plan} / M_{par}$')
ax[axNum].set_xlim(xLim1, doPlan1.tMax)
ax[axNum].axvline(t_vline1, color=(0,0,1,1), linestyle='--')
ax[axNum].axvline(t_vline2, color=(1,0,0,1), linestyle='--')

# plot nClumps over time
axNum = 2
ax[axNum].plot(doPlan1.time[nStart:], doPlan1.nClumpsList[nStart:], 'k', linewidth=2)
ax[axNum].set_ylabel(r'$N_{clumps}$')
ax[axNum].set_xlabel(r'$t \Omega$')
ax[axNum].set_xlim(xLim1, doPlan1.tMax)
ax[axNum].axvline(t_vline1, color=(0,0,1,1), linestyle='--')
ax[axNum].axvline(t_vline2, color=(1,0,0,1), linestyle='--')

# close and save figure
plt.tight_layout()
plt.savefig(pathSave + "p_vs_t.png", bbox_inches='tight')
plt.close('all')

################################################################################












#
