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
import corner
################################################################################
# paths and CL args
pathBase = str(sys.argv[1])
nStart   = int(sys.argv[2])
nEnd     = int(sys.argv[3])
pathPlan = pathBase + 'planOutput2/'
pathSave = pathBase + 'plots/plan_mcmcFits/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
doPlan = readerPlan.DataPlan(pathPlan, nStart=nStart, nTot=nEnd, nPar=512*512*512)
################################################################################

spl_master  = np.zeros([nEnd-nStart+1, 3, 1])
stpl_master = np.zeros([nEnd-nStart+1, 3, 2])
tpl_master  = np.zeros([nEnd-nStart+1, 3, 2])
bcpl_master = np.zeros([nEnd-nStart+1, 3, 3])
bpl_master  = np.zeros([nEnd-nStart+1, 3, 3])
vtpl_master = np.zeros([nEnd-nStart+1, 3, 3])
tspl_master = np.zeros([nEnd-nStart+1, 3, 5])

for n in range(nStart, nEnd):
    print("\n\n\n\ndoing all fits for n=" + str(n))
    mp1, ngtm = readerPlan.getCumMassHist(doPlan, n)
    # SPL
    spl_fit = pstats.Fit_Pipeline(mp1, pstats.fitInfo_spl)
    spl_fit.cp2.savefig(pathSave + "cp_"+ pstats.fitInfo_spl.name + "_" + str(n) + ".png")
    # STPL
    stpl_fit = pstats.Fit_Pipeline(mp1, pstats.fitInfo_stpl)
    stpl_fit.cp2.savefig(pathSave + "cp_"+ pstats.fitInfo_stpl.name + "_" + str(n) +".png")
    # TPL
    tpl_fit = pstats.Fit_Pipeline(mp1, pstats.fitInfo_tpl)
    tpl_fit.cp2.savefig(pathSave + "cp_"+ pstats.fitInfo_tpl.name + "_" + str(n) +".png")
    # BCPL
    bcpl_fit = pstats.Fit_Pipeline(mp1, pstats.fitInfo_bcpl)
    bcpl_fit.cp2.savefig(pathSave + "cp_"+ pstats.fitInfo_bcpl.name + "_" + str(n) +".png")
    # BPL
    bpl_fit = pstats.Fit_Pipeline(mp1, pstats.fitInfo_bpl)
    bpl_fit.cp2.savefig(pathSave + "cp_"+ pstats.fitInfo_bpl.name + "_" + str(n) +".png")
    # VTPL
    vtpl_fit = pstats.Fit_Pipeline(mp1, pstats.fitInfo_vtpl)
    vtpl_fit.cp2.savefig(pathSave + "cp_"+ pstats.fitInfo_vtpl.name + "_" + str(n) +".png")
    # TSPL
    tspl_fit = pstats.Fit_Pipeline(mp1, pstats.fitInfo_tspl)
    tspl_fit.cp2.savefig(pathSave + "cp_"+ pstats.fitInfo_tspl.name + "_" + str(n) +".png")

    plt.clf()

    # assign fit params to master arrays
    nEff = n-nStart
    spl_master[nEff]  = spl_fit.params_mcmc
    stpl_master[nEff] = stpl_fit.params_mcmc
    tpl_master[nEff]  = tpl_fit.params_mcmc
    bcpl_master[nEff] = bcpl_fit.params_mcmc
    bpl_master[nEff]  = bpl_fit.params_mcmc
    vtpl_master[nEff] = vtpl_fit.params_mcmc
    tspl_master[nEff] = tspl_fit.params_mcmc

np.save(pathSave + "params_spl.npy", spl_master)
np.save(pathSave + "params_stpl.npy", stpl_master)
np.save(pathSave + "params_tpl.npy", tpl_master)
np.save(pathSave + "params_bcpl.npy", bcpl_master)
np.save(pathSave + "params_bpl.npy", bpl_master)
np.save(pathSave + "params_vtpl.npy", vtpl_master)
np.save(pathSave + "params_tspl.npy", tspl_master)



#
