#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import athenaReader1d as reader1d
import athenaReaderPhst as readerPhst
import athenaTools as tools
from matplotlib.backends.backend_pdf import PdfPages
################################################################################
if str(sys.argv[1])=='10':
	pathBase = '../../data/prodRuns/'
	runNameList = ['run10','run11','run12','run13']
	alphaInList = [1.e-3, 1.e-4, 1.e-5, 0.0]
	tsList      = [0.3, 0.3, 0.3, 0.3]
	colorList   = ['tab:blue', 'tab:orange', 'g', 'r']
	pathSave = pathBase + 'plots/turbAnalysis10/'
################################################################################
elif str(sys.argv[1])=='2?':
	pathBase = '../../data/prodRuns/'
	runNameList = ['run20','run21','run22','run23']
	alphaInList = [1.e-3, 1.e-4, 1.e-5, 9.e-6]
	tsList      = [0.3, 0.3, 0.3, 0.3]
	colorList   = ['tab:blue', 'tab:orange', 'g', 'r']
	pathSave = pathBase + 'plots/turbAnalysis20/'
################################################################################
elif str(sys.argv[1])=='mix':
	pathBase = '../../data/prodRuns/'
	runNameList = ['run10','run20','run30']
	alphaInList = [1.e-3, 1.e-3, 1.e-3]
	tsList      = [0.3, 0.3, 0.3, 0.3]
	colorList   = ['tab:blue', 'tab:orange', 'g', 'r']
	pathSave = pathBase + 'plots/turbAnalysisMix/'
################################################################################
elif str(sys.argv[1])=='6?':
	pathBase = '../../data/prodRuns/'
	runNameList = ['run60','run61','run62','run63']
	alphaInList = [1.e-3, 1.e-4, 1.e-5, 0.0]
	tsList      = [0.3, 0.3, 0.3, 0.3]
	colorList   = ['tab:blue', 'tab:orange', 'g', 'r']
	pathSave = pathBase + 'plots/turbAnalysis60/'
################################################################################
elif str(sys.argv[1])=='7?':
	pathBase = '../../data/prodRuns/'
	runNameList = ['run70','run71','run72']
	alphaInList = [1.e-3, 1.e-3, 1.e-3]
	tsList      = [0.3, 0.3, 0.3]
	colorList   = ['tab:blue', 'tab:orange', 'g', 'r']
	pathSave = pathBase + 'plots/turbAnalysis70/'
################################################################################
elif str(sys.argv[1])=='8?':
	pathBase = '../../data/prodRuns/'
	runNameList = ['run80','run81','run82']
	alphaInList = [1.e-3, 1.e-3, 1.e-3]
	tsList      = [0.3, 0.3, 0.3]
	colorList   = ['tab:blue', 'tab:orange', 'g', 'r']
	pathSave = pathBase + 'plots/turbAnalysis80/'
################################################################################
elif str(sys.argv[1])=='9?':
	pathBase = '../../data/prodRuns/'
	runNameList = ['run90','run91']#,'run92']
	alphaInList = [1.e-3, 1.e-3]#, 1.e-3]
	tsList      = [0.3, 0.3]#, 0.3]
	colorList   = ['tab:blue', 'tab:orange', 'g', 'r']
	pathSave = pathBase + 'plots/turbAnalysis90/'
################################################################################
elif str(sys.argv[1])=='11?':
	pathBase = '../../data/prodRuns/'
	runNameList = ['run110','run111','run112']
	alphaInList = [1.e-3, 1.e-3, 1.e-3]
	tsList      = [0.3, 0.3, 0.3]
	colorList   = ['tab:blue', 'tab:orange', 'g', 'r']
	pathSave = pathBase + 'plots/turbAnalysis110/'
################################################################################
elif str(sys.argv[1])=='7and8':
	pathBase = '../../data/prodRuns/'
	runNameList = ['run70', 'run80', 'run71', 'run81', 'run72', 'run82']
	alphaInList = [1.e-3, 1.e-3, 1.e-3, 1.e-3, 1.e-3, 1.e-3]
	tsList      = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
	a           = 0.4
	colorList   = [(1,0,0,1), (1,0,0,a), (0.5,1,0,1), (0.5,1,0,a), (0,0,1,1), (0,0,1,a)]
	pathSave = pathBase + 'plots/turbAnalysis7and8/'
################################################################################

if not os.path.exists(pathSave): os.makedirs(pathSave)

# read 1d and phst files in for all runs
do1dList    = []
doPhstList  = []
for n in range(len(runNameList)):
	path1d     = pathBase + runNameList[n] + '/1d/'
	pathParHst = pathBase + runNameList[n] + '/'
	do1dList  .append(reader1d.Data1d(path1d))
	doPhstList.append(readerPhst.DataPhst(pathParHst))
	do1dList[n].addCol(reader1d.dv,     'dv',     r'$\delta v$')
	do1dList[n].addCol(reader1d.dvz,    'dvz',    r'$\delta v_z$')
	do1dList[n].addCol(reader1d.alpha,  'alpha',  r'$\alpha_z$')
	do1dList[n].addCol(reader1d.alphaz, 'alphaz', r'$\alpha$')




# plot dv over time vs expected value based on alpha_in
key = 'dv'
for n in range(len(do1dList)):
	do        = do1dList[n]
	alphaIn   = alphaInList[n]
	color     = colorList[n]
	reader1d.timeEvo(do, key, legendLabel=r'$\alpha_{in}=$'+str(alphaIn), logForce=1, color=color)
	plt.axhline(y=np.sqrt(alphaIn), linestyle='--', color=color)
plt.legend()
tools.saveAndClear(pathSave + "gas_timeEvo_" + key + ".png")

# plot dv at midplane over time vs expected value based on alpha_in
key = 'dv'
for n in range(len(do1dList)):
	do        = do1dList[n]
	alphaIn   = alphaInList[n]
	color     = colorList[n]
	reader1d.timeEvo(do, key, z1=-do.zmax/4.0, z2=do.zmax/4.0,
					 legendLabel=r'$\alpha_{in}=$'+str(alphaIn), logForce=1, color=color)
	plt.axhline(y=np.sqrt(alphaIn), linestyle='--', color=color)
plt.legend()
tools.saveAndClear(pathSave + "gas_timeEvo_mid_" + key + ".png")





# plot dvz over time vs expected value based on alpha_in
key = 'dvz'
for n in range(len(do1dList)):
	do        = do1dList[n]
	alphaIn   = alphaInList[n]
	color     = colorList[n]
	reader1d.timeEvo(do, key, legendLabel=r'$\alpha_{in}=$'+str(alphaIn), logForce=1, color=color)
	plt.axhline(y=np.sqrt(alphaIn), linestyle='--', color=color)
plt.legend()
tools.saveAndClear(pathSave + "gas_timeEvo_" + key + ".png")

# plot dvz at midplane over time vs expected value based on alpha_in
key = 'dvz'
for n in range(len(do1dList)):
	do        = do1dList[n]
	alphaIn   = alphaInList[n]
	color     = colorList[n]
	reader1d.timeEvo(do, key, z1=-do.zmax/4.0, z2=do.zmax/4.0,
					 legendLabel=r'$\alpha_{in}=$'+str(alphaIn), logForce=1, color=color)
	plt.axhline(y=np.sqrt(alphaIn), linestyle='--', color=color)
plt.legend()
tools.saveAndClear(pathSave + "gas_timeEvo_mid_" + key + ".png")









# plot particle dz vs expected value based on alpha ~ dvz^2
for n in range(len(do1dList)):
	do         = do1dList[n]
	alphaIn    = alphaInList[n]
	color      = colorList[n]
	doPhst     = doPhstList[n]
	readerPhst.timeEvo(doPhst, 'zvar', legendLabel=r'$\alpha_{in}=$'+str(alphaIn), logForce=1)
	alphazMean = np.mean(do.data['alphaz'][do.nt//2:])
	parh       = np.sqrt(alphazMean/tsList[n])
	#plt.axhline(y=parh, linestyle='--', color=color)
#plt.axhline(y=do.dz, linestyle='--', color='k')
#plt.axhline(y=do.zmax*2.0, linestyle='--', color='k')
plt.legend()
tools.saveAndClear(pathSave + "par_" + 'scaleHeightComparison' + ".png")

# plot particle dz vs expected value based on alpha ~ dvz^2, just at mid-plane
for n in range(len(do1dList)):
	do         = do1dList[n]
	alphaIn    = alphaInList[n]
	color      = colorList[n]
	doPhst     = doPhstList[n]
	alphazMean = np.mean(do.data['alphaz'][do.nt//2:, do.nz//2-8:do.nz//2+8])
	parh       = np.sqrt(alphazMean/tsList[n])
	readerPhst.timeEvo(doPhst, 'zvar', legendLabel=r'$\alpha_{in}=$'+str(alphaIn), logForce=1)
	#plt.axhline(y=parh, linestyle='--', color=color)
	if alphaIn == 1.e-3:
		print(alphaIn)
		print(np.round(alphazMean,6))
		print(np.round(parh,6))
		print(np.round(np.mean(doPhst.data['zvar'][doPhst.nt//2:]),4))
#plt.axhline(y=do.dz, linestyle='--', color='k')
#plt.axhline(y=do.zmax*2.0, linestyle='--', color='k')
plt.legend()
tools.saveAndClear(pathSave + "par_" + 'scaleHeightComparison2' + ".png")

# plot measured alpha vs input alpha
alphaActList  = []
for n in range(len(do1dList)):
	do      = do1dList[n]
	alphaActList.append(np.mean(do.data['alpha'][do.nt//2:,do.nz//2-10:do.nz//2+10]))
plt.loglog(alphaInList, alphaActList, 'ko', label='data')
plt.xlabel(r'$\alpha_{in}$')
plt.ylabel(r'$\alpha_{real}$')
alphaSat = [1.e-8, 1.e-7, 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1.e0, 1.e1]
plt.loglog(alphaSat, alphaSat, lineStyle='--', color='r', label='dissipation')
alphaIn  = []
for x in alphaSat:
	alphaIn.append(np.power(100.0*x, 1.5))
plt.loglog(alphaIn, alphaSat, lineStyle='--', color='g', label='outflow')
BoA      = 95
alphaIn  = []
for x in alphaSat:
	alphaIn.append(x + BoA*np.power(x, 1.5))
plt.loglog(alphaIn, alphaSat, lineStyle='--', color='b', label='both')
plt.legend()
plt.xlim(8.e-9, 2.e1)
tools.saveAndClear(pathSave + "alphaAct.png")
