#!/usr/bin/python
import numpy as np
import matplotlib as m
#m.use('Agg')
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import plan_stats as pstats
import athenaReader1d as reader1d
import athenaReaderPhst as readerPhst
import athenaTools as tools
import planOutputReader as readerPlan
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import numpy.polynomial.polynomial as poly
################################################################################
# paths and CL args
pathBase    = str(sys.argv[1])
pathFitData = pathBase + 'plots/plan_mcmcFits/'
pathSave    = pathBase + 'plots/plan_mcmcFitPlots/'
pathPlan    = pathBase + 'planOutput2/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
tStart = 23.0

fitInfoList = [pstats.fitInfo_spl,
               pstats.fitInfo_tpl,
               pstats.fitInfo_stpl,
               pstats.fitInfo_vtpl,
               pstats.fitInfo_bcpl,
               pstats.fitInfo_bpl,
               pstats.fitInfo_tspl]

for fitInfo in fitInfoList:
    params = np.load(pathFitData + "params_" + fitInfo.name + ".npy")
    nPlots = params.shape[2]
    nt     = params.shape[0]
    timeArr = np.arange(nt)/10.0 + tStart

    fig = plt.figure(figsize=(10, 3*nPlots), dpi=120)
    ax  = []
    for n in range(nPlots):
        ax.append(plt.subplot2grid((nPlots, 1), (n, 0), rowspan=1, colspan=1))
    for n in range(nPlots):
        ax[n].plot(timeArr, params[:,0,n], 'k')
        ax[n].fill_between(timeArr, params[:,1,n], params[:,2,n], color='gray', alpha=0.5)
        ax[n].set_ylabel(fitInfo.paramNames[n])
        ax[n].set_xlabel(r'$t \Omega$')
    plt.savefig(pathSave + fitInfo.name + ".png", bbox_inches='tight'); plt.clf();




























#
