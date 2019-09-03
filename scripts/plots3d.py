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
alpha    = float(sys.argv[4])
dvBase   = np.sqrt(alpha)
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/plots3d/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
plt.figure(0)
################################################################################
key = 'dpar'
for n in [5,15,25,35]:
	print(n)
	reader3d.profile(do3d, key, color='r', figNum=0, nStart=n, nEnd=n+5, legendLabel=do3d.header[key])
	plt.ylim(1.e-3, 1.e1)
	mean = 0
	for i in range(5):
		print(n+i)
		data   =   do3d.get3d(key, n+i)
		if i==0: data1d  = np.mean(data, axis=(0,1))
		else:    data1d += np.mean(data, axis=(0,1))
	data1d /= 5
	plt.title(np.amax(data)) 
	tools.saveAndClear(pathSave + 'profile_dpar_'+str(n)+'_'+str(n+5)+'.png', figNum=0)


# dv profiles, all components on same plot
key = 'dv'
reader3d.profile(do3d, key, color='k', figNum=0, nStart=nStart, nEnd=nEnd, legendLabel=do3d.header[key])
tools.saveAndClear(pathSave + 'profilePert_newdv.png', figNum=0)
key = 'dvx'
reader3d.profile(do3d, key, color='r', figNum=0, nStart=nStart, nEnd=nEnd, legendLabel=do3d.header[key])
key = 'dvy'
reader3d.profile(do3d, key, color='g', figNum=0, nStart=nStart, nEnd=nEnd, legendLabel=do3d.header[key])
key = 'dvz'
reader3d.profile(do3d, key, color='b', figNum=0, nStart=nStart, nEnd=nEnd, legendLabel=do3d.header[key])
plt.ylim(dvBase*0.1, dvBase*3.2)
plt.axhline(dvBase, color='gray', linestyle='--')
#plt.axhline(dvBase/np.sqrt(3.), color='gray', linestyle=':')
plt.legend(loc='best')
tools.saveAndClear(pathSave + 'profilePert_newdv.png', figNum=0)

################################################################################
# dv time evo, all components on same plot
'''
key = 'dv'
reader3d.timeEvo(do3d, key, color='k', figNum=0, legendLabel=do3d.header[key])
tools.saveAndClear(pathSave + 'timeEvoPert_newdv.png', figNum=0)
key = 'dvx'
reader3d.timeEvo(do3d, key, color='r', figNum=0, legendLabel=do3d.header[key])
key = 'dvy'
reader3d.timeEvo(do3d, key, color='g', figNum=0, legendLabel=do3d.header[key])
key = 'dvz'
reader3d.timeEvo(do3d, key, color='b', figNum=0, legendLabel=do3d.header[key])
plt.ylim(dvBase*0.1, dvBase*3.2)
plt.axhline(dvBase, color='gray', linestyle='--')
#plt.axhline(dvBase/np.sqrt(3.), color='gray', linestyle=':')
plt.legend(loc='best')
tools.saveAndClear(pathSave + 'timeEvoPert_alldv.png', figNum=0)
'''
################################################################################
# max dpar

plt.figure(0)
print("making max dpar plot...")
sys.stdout.flush()
key = 'dpar'
plotDataList = []
for n in range(do3d.nt):
	print("reading in time step " + str(n) + "...")
	sys.stdout.flush()
	data = do3d.get3d(key, n)
	plotDataList.append(np.amax(data))
plotData = np.asarray(plotDataList)
print(plotData)
plt.xlabel(r'$t \Omega$')
plt.ylabel('MAX(' + do3d.header[key] + ')')
plt.semilogy(do3d.t, plotData, color='k')
plt.axvline(20, color='gray', lineStyle=':')
plt.axhline(plotData[20], color='gray', lineStyle=':')
plt.xlim(0.0, np.amax(do3d.t))
tools.saveAndClear(pathSave + 'timeEvo_maxDPar.png', figNum=0)

################################################################################
'''
# dv distribution
key = 'dv'
n   = nEnd
data3d = do3d.get3d(key, n)
i1=-5; i2=0;
bins = np.logspace(i1,i2,num=200)
plt.hist(data3d.flatten(), bins=bins, log=True, color=(0,0,0,0.3), density=True)
plt.xscale('log')
plt.xlim(10**i1,10**i2)
#plt.ylim(1.0, 1.e4*do3d.nx*do3d.ny*do3d.nz/(64.*64.*64.))
plt.ylabel('Probability Density')
plt.xlabel(do3d.header[key])
tools.saveAndClear(pathSave + 'timeEvo_maxDPar.png', figNum=0)
'''






















# perts
#for key in ['dv', 'dvx', 'dvy', 'dvz', 'drho']:
    #reader3d.profile(do3d, key, figNum=0, nStart=nStart, nEnd=nEnd)
    #plt.ylim(3.e-3, 1.e-1)
    #tools.saveAndClear(pathSave + 'profilePert_' + key + '.png', figNum=0)
    #reader3d.timeEvo(do3d, key, figNum=0)
    #tools.saveAndClear(pathSave + 'timeEvoPert_' + key + '.png', figNum=0)
################################################################################

# simple averages of quantities
#for key in ['vx', 'vy', 'vz']:
    #reader3d.profile(do3d, key, figNum=0, absAvg=0, absPlot=0)
    #tools.saveAndClear(pathSave + 'profileRealAvg_' + key + '.png', figNum=0)
    #reader3d.timeEvo(do3d, key, figNum=0, absAvg=0, absPlot=0)
    #tools.saveAndClear(pathSave + 'timvoRealAvg_' + key + '.png', figNum=0)
################################################################################
# abs averages of quantities
#for key in ['vx', 'vy', 'vz']:
    #reader3d.profile(do3d, key, figNum=0)
    #tools.saveAndClear(pathSave + 'profileAbsAvg_' + key + '.png', figNum=0)
    #reader3d.timeEvo(do3d, key, figNum=0)
    #tools.saveAndClear(pathSave + 'timeEvoAbsAvg_' + key + '.png', figNum=0)
################################################################################
# normalized perts
#for key in ['drhoNorm', 'dvxNorm', 'dvyNorm', 'dvzNorm', 'dvNorm']:
    #reader3d.profile(do3d, key, figNum=0)
    #tools.saveAndClear(pathSave + 'profilePertNorm_' + key + '.png', figNum=0)
    #reader3d.timeEvo(do3d, key, figNum=0)
    #tools.saveAndClear(pathSave + 'timeEvoPertNorm_' + key + '.png', figNum=0)
################################################################################

















#
