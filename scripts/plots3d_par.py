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
import multiprocessing as mp
################################################################################
pathBase = str(sys.argv[1])
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/plots3d/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
################################################################################
# simple averages of quantities
def realAvgProfile(do3d, key, pathSave):
	reader3d.profile(do3d, key, figNum=0, absAvg=0, absPlot=0)
	tools.saveAndClear(pathSave + 'profileRealAvg_' + key + '.png', figNum=0)
def realAvgTimeEvo(do3d, key, pathSave):
	reader3d.timeEvo(do3d, key, figNum=0, absAvg=0, absPlot=0)
	tools.saveAndClear(pathSave + 'timeEvoRealAvg_' + key + '.png', figNum=0)
def absAvgProfile(do3d, key, pathSave):
	reader3d.profile(do3d, key, figNum=0)
	tools.saveAndClear(pathSave + 'profileAbsAvg_' + key + '.png', figNum=0)
def absAvgTimeEvo(do3d, key, pathSave):
	reader3d.timeEvo(do3d, key, figNum=0)
	tools.saveAndClear(pathSave + 'timeEvoAbsAvg_' + key + '.png', figNum=0)
def pertProfile(do3d, key, pathSave):
	reader3d.profile(do3d, key, figNum=0)
	tools.saveAndClear(pathSave + 'profilePert_' + key + '.png', figNum=0)
def pertTimeEvo(do3d, key, pathSave):
	reader3d.timeEvo(do3d, key, figNum=0)
	tools.saveAndClear(pathSave + 'timeEvoPert_' + key + '.png', figNum=0)
################################################################################

jobList = []
for key in ['vx', 'vy', 'vz']:
	job = mp.Process(target=realAvgProfile, args=(do3d, key, pathSave))
	jobList.append(job)
	job = mp.Process(target=realAvgTimeEvo, args=(do3d, key, pathSave))
	jobList.append(job)
	job = mp.Process(target=absAvgProfile, args=(do3d, key, pathSave))
	jobList.append(job)
	job = mp.Process(target=absAvgTimeEvo, args=(do3d, key, pathSave))
	jobList.append(job)

for key in ['drho', 'dvx', 'dvy', 'dvz']:
	job = mp.Process(target=pertProfile, args=(do3d, key, pathSave))
	jobList.append(job)
	job = mp.Process(target=pertTimeEvo, args=(do3d, key, pathSave))
	jobList.append(job)

for job in jobList:
	job.start()


reader3d.profile(do3d, 'dvz', figNum=0, absAvg=1, absPlot=1)
tools.saveAndClear(pathSave + 'testing1.png', figNum=0)
plotData = np.zeros(do3d.nt)
for n in range(len(plotData)):
	data = do3d.get3d('dvz', n)
	plotData[n] = np.mean(np.absolute(data))
plt.plot(do3d.t, plotData)
tools.saveAndClear(pathSave + 'testing2.png', figNum=0)












################################################################################
# normalized perts
#for key in ['drhoNorm', 'dvxNorm', 'dvyNorm', 'dvzNorm', 'dvNorm']:
    #reader3d.profile(do3d, key, figNum=0)
    #tools.saveAndClear(pathSave + 'profilePertNorm_' + key + '.png', figNum=0)
    #reader3d.timeEvo(do3d, key, figNum=0)
    #tools.saveAndClear(pathSave + 'timeEvoPertNorm_' + key + '.png', figNum=0)
################################################################################













#
