#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import math
import sys
import time
sys.path.append('../python')
import athenaReader3d as reader3d
import athenaTools as tools
from matplotlib.backends.backend_pdf import PdfPages
import multiprocessing as mp
################################################################################
pathBase = str(sys.argv[1])
myNpc    = int(sys.argv[2])
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/slices/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
plt.figure(0)
################################################################################
def oneSlice(do3d, key, n):
    plt.figure(0)
    t0 = time.time()
    reader3d.slicePlot(do3d, key, n=n, figNum=0)
    tools.saveAndClear(pathSave + 'midplaneSlice_' + key + '_' + str(n) + '.png', figNum=0)
    t = np.round(time.time() - t0, 2)
    print('this plot took ' + str(t) + 's')
    sys.stdout.flush()
################################################################################
jobList = []
for key in ['drho','dpar','dvx','dvy','dvz']:
    for n in range(0, do3d.nt, 1):
        job = mp.Process(target=oneSlice, args=(do3d, key, n))
        jobList.append(job)
while len(jobList)>0:
    theseJobs = []
    for n in range(myNpc):
        try:    theseJobs.append(jobList.pop(0))
        except: a=1
    print('################################################################')
    t0 = time.time()
    print('starting up ' + str(myNpc) + ' processes')
    sys.stdout.flush()
    for job in theseJobs:
        job.start()
    for job in theseJobs:
        job.join()
    t1 = time.time() - t0
    print('these ' + str(myNpc) + ' processes are ending after ' + str(t1))
    print('################################################################\n')
    sys.stdout.flush()





















#
