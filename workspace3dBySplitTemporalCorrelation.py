#!/usr/bin/python
from __future__ import unicode_literals
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import athenaReader1d as reader1d
import athenaReader3d as reader3d
import os
from matplotlib.backends.backend_pdf import PdfPages
import resource
m.rcParams['text.usetex'] = True
m.rcParams['text.latex.unicode'] = True
font = {'family' : 'normal', 'weight' : 'bold', 'size'   : 14}
m.rc('font', **font)

def getTimeStepString(i):
	if i > 999:	 zstring = ""
	elif i > 99: zstring = "0"
	elif i > 9:  zstring = "00"
	elif i > -1: zstring = "000"
	return zstring+str(i)

upScaleFactor = 10.0
def getCorrArray(aStatic, aShift, shiftMax):
	# try a range of shifts and calculate difference in curves
	corrArray1  = np.zeros(2*shiftMax*int(upScaleFactor))
	shiftArray  = np.arange(-shiftMax, shiftMax, 1.0)
	shiftArray1 = np.arange(-shiftMax, shiftMax, 1.0/upScaleFactor)
	shiftIndexArray1 = np.arange(int(-shiftMax*upScaleFactor), int(shiftMax*upScaleFactor), 1)
	timeArray   = np.arange(0, len(aStatic), 1.0)
	timeArray1  = np.arange(0, len(aStatic), 1.0/upScaleFactor)
	aStatic1    = np.interp(timeArray1, timeArray, aStatic)   		
	aShift1     = np.interp(timeArray1, timeArray, aShift)
	nShift      = 0
	for shiftIndex in shiftIndexArray1:
		# define comparison array and shift shifted array		
		unshifted = aStatic1
		shifted   = np.roll(aShift1, shiftIndex)
		# crop to get rid of edge effects
		unshifted = unshifted[int(shiftMax*upScaleFactor):int(-shiftMax*upScaleFactor)]
		shifted   = shifted  [int(shiftMax*upScaleFactor):int(-shiftMax*upScaleFactor)]
		corr 			= np.sum(shifted*unshifted)
		corrArray1[nShift] = corr
		nShift+=1
	return shiftArray1, corrArray1

def getShiftMaxCoord(shiftArray, a):
	arg = np.argmax(a)
	#corrNorm = a[arg] / np.mean(np.absolute(a))
	corrNorm = a[arg]
	return corrNorm, shiftArray[arg]
	#return shiftArray[300+np.argmax(a[300:700])]

amDict = {'4100'     : 1.e-2,
          '4101'     : 1.8e-2,
          '4102'     : 3.2e-2,
          '4103'     : 5.6e-2,
          '4104'     : 1.e-1,
          '4105'     : 1.8e-1,
          '4106'     : 3.2e-1,
          '4107'     : 5.6e-1,
          '4108'     : 1.e0,
          '4109'     : 1.8e0,
          '4110'     : 3.2e0,
          '4111'     : 5.6e0,
          '4112'     : 1.e1  }

midPlaneLagList = np.zeros(13)
coronaLagList = np.zeros(13)
betaMidList     = np.zeros(13)
betaMinMidList  = np.zeros(13)
betaTransList   = np.zeros(13)
betaCoronaList  = np.zeros(13)
midCorrList  = np.zeros(13)
betaProfileList = np.zeros([13,30*8])
nRun = 1

nSlices = 240
for idIn in sys.argv[1:]:
	nRun +=1
	idString = str(idIn)
	idInt    = int(idIn)
	print(idString)

	suffix = '.By.npy'
	dt = 1.0
	path = '/data/adMedResRuns/3d/combTabsSplit'+idString+'/'
	pdf = PdfPages(path+"by_timeCorrHighRes.pdf")
	#pdf = PdfPages(path+"by_timeCorr.pdf")
	baseName = 'StratCooling'

	tTemp = np.arange(0, 10000, dt)
	x  = np.load(path+baseName+'.0000.x.npy'); nx = x.shape[0];
	y  = np.load(path+baseName+'.0000.y.npy'); ny = y.shape[0]; 
	z  = np.load(path+baseName+'.0000.z.npy'); nz = z.shape[0]; 
	dx = x[1]-x[0]

	ntMax = 350
	dataTemp = np.zeros([1000, x.shape[0], y.shape[0], z.shape[0]])
	for n in range(1000):
	#for n in range(200):
		fileName = path+baseName+'.'+getTimeStepString(n)+'.By.npy'
		print('reading in ' + fileName)
		if os.path.exists(fileName) and n!=ntMax:
			loadedData = np.load(fileName)
			dataTemp[n] = loadedData
		else:
			nt = n
			break
	#nt=100

	data = dataTemp[0:nt,:,:,:]
	t = tTemp[0:nt]
	
	# ST plot
	aspect   = 0.2*(t[-1])/z[-1]
	extent   = [0, t[-1], -z[-1], z[-1]]
	plotData = np.mean(data, axis=(1,2))
	plotData = np.transpose(np.fliplr(plotData))
	maxAbsVal = np.amax(np.absolute(plotData))
	plt.imshow(plotData, extent=extent, aspect=aspect, cmap='coolwarm', 
						 norm=colors.SymLogNorm(linthresh=maxAbsVal/100, linscale=0.5, 
																		vmin=-maxAbsVal, vmax=maxAbsVal))
	plt.title(r"$B_y$"); plt.xlabel("time (orbits)"); plt.ylabel(r"$z/H$");
	plt.colorbar(shrink=0.5); plt.tight_layout();
	plt.savefig(pdf, format='pdf', bbox_inches='tight');
	plt.savefig(path+"st_by.png", bbox_inches='tight'); plt.clf();

	# spatial averages vs time
	plotData1 = np.mean(data[:,:,:,int(3.2/dx):int(4.8/dx)], axis=(1,2,3))
	plotData2 = np.mean(data[:,:,:,int(1.5/dx):int(2.5/dx)], axis=(1,2,3))
	plotData3 = np.mean(data[:,:,:,int(5.5/dx):int(6.5/dx)], axis=(1,2,3))
	plt.plot(t, plotData1, label='midplane')
	plt.plot(t, plotData2, label='-corona')
	plt.plot(t, plotData3, label='+corona')
	plt.xlabel("time (orbits)"); plt.ylabel(r"$<B_y>$");
	plt.legend(loc=(1.02, 0.0))
	plt.savefig(pdf, format='pdf', bbox_inches='tight');
	plt.savefig(path+"byVsTime.png", bbox_inches='tight'); plt.clf();

	zTrans = 1.12
	ziTrans1     = int((4-zTrans)/dx)
	ziTrans2     = int((4+zTrans)/dx) + 1

	#zRef = 0.95
	zRef = 0.0
	ziRef1     = int((4-zRef)/dx)
	ziRef2     = int((4+zRef)/dx) + 1

	slices = [];
	step = nz/nSlices;
	for i in range(nSlices):
		zi1  = int(i*step)
		zi2  = int((i+1)*step)
		slices.append(np.mean(data[:,:,:,zi1:zi2], axis=(1,2,3)))

	slicesRand = [];
	for i in range(nSlices):
		zi1  = int(i*step)
		zi2  = int((i+1)*step)
		slicesRand.append((np.random.rand(nt)*2.0-1.0)*np.mean(np.absolute(slices[i])))
		#plt.plot(slicesRand[i]); plt.show(); plt.clf();

	lags      = np.zeros(nSlices)
	corrsMax  = np.zeros(nSlices)
	corrsRand = np.zeros(nSlices)
	allCorrs  = np.zeros([nSlices,400])
	# compute correlations

	# +side correlations with first +side slice
	for i in range(nSlices//2,nSlices):
		shiftArray, corrArray = getCorrArray(slices[ziRef2], slices[i], 20)
		allCorrs[i] = corrArray
		if np.mod(i,5)==0:
			plt.plot(shiftArray, corrArray, label = str(nSlices//2) + ", " + str(i));
		corrsMax[i], lags[i] = getShiftMaxCoord(shiftArray, corrArray)
		shiftArray, corrArray = getCorrArray(slices[nSlices//2], slicesRand[i], 20)
		corrsRand[i] = np.amax(corrArray)
	plt.axhline(y=0, linestyle='--', color='k')
	plt.axvline(x=0, linestyle='--', color='k')
	plt.xlabel('shift'); plt.ylabel('correlation'); plt.legend(loc=(1.02,0.0));
	plt.savefig(pdf, format='pdf', bbox_inches='tight');
	plt.clf();

	# -side correlations with first -side slice
	for ii in range(0,nSlices//2):
		i = nSlices//2-1-ii
		shiftArray, corrArray = getCorrArray(slices[ziRef1], slices[i], 20)
		allCorrs[ii] = corrArray
		if np.mod(ii,5)==0:
			plt.plot(shiftArray, corrArray, label = str(nSlices//2-1) + ", " + str(i));
		corrsMax[i], lags[i] = getShiftMaxCoord(shiftArray, corrArray)
		shiftArray, corrArray = getCorrArray(slices[nSlices//2-1], slicesRand[i], 20)
		corrsRand[i] = np.amax(corrArray)
	plt.axhline(y=0, linestyle='--', color='k')
	plt.axvline(x=0, linestyle='--', color='k')
	plt.xlabel('shift'); plt.ylabel('correlation'); plt.legend(loc=(1.02,0.0));
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();

	midCorrList[nRun] = np.median(corrsMax[ziTrans1:ziTrans2]) / np.amax(corrsMax)
	
	# ST plot of correlations
	aspect   = 0.2*(t[-1])/z[-1]
	extent   = [0, t[-1], -z[-1], z[-1]]
	plotData = allCorrs
	#plotData = np.transpose(np.fliplr(plotData))
	maxAbsVal = np.amax(np.absolute(plotData))
	plt.imshow(plotData, extent=extent, aspect=aspect, cmap='coolwarm', 
						 norm=colors.SymLogNorm(linthresh=maxAbsVal/100, linscale=0.5, 
																		vmin=-maxAbsVal, vmax=maxAbsVal))
	plt.title(r"$B_y$"); plt.xlabel("time (orbits)"); plt.ylabel(r"$z/H$");
	plt.colorbar(shrink=0.5); plt.tight_layout();
	plt.savefig(pdf, format='pdf', bbox_inches='tight');
	plt.savefig(path+"st_by.png", bbox_inches='tight'); plt.clf();


	plt.plot(z, lags, 'bo')
	plt.xlabel('z/H');
	plt.axvline(x=-zTrans, linestyle='--', color='k'); plt.axvline(x=zTrans, linestyle='--', color='k');
	plt.axvline(x=-zRef, linestyle='--', color='k');   plt.axvline(x=zRef, linestyle='--', color='k');
	midPlaneLag1 = np.amax(lags[nz//2-45:nz//2]) 
	midPlaneLag2 = np.amax(lags[nz//2:nz//2+45]) 
	midPlaneLag = 0.5 * (midPlaneLag1 + midPlaneLag2)
	midPlaneLagList[nRun] = midPlaneLag
	coronaLag1 = np.mean(lags[0*30:1*30]) 
	coronaLag2 = np.mean(lags[7*30:8*30]) 
	coronaLag = 0.5 * (coronaLag1 + coronaLag2)
	coronaLagList[nRun] = np.absolute(coronaLag)
	print(coronaLagList)
	plt.ylabel('most correlated time shift (orbits)');
	plt.axhline(y=0, linestyle='--', color='k') 
	plt.savefig(path+"lagVsZ.png", bbox_inches='tight');
	plt.title("Am=" + str(amDict[idString]) + ', ' + str(midPlaneLag))
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();


	corrsNorm = corrsMax / corrsRand
	binFactor = 1
	binCorrsNorm = np.zeros(nz//binFactor)
	for i in range(len(binCorrsNorm)):
		binCorrsNorm[i] = np.mean(corrsNorm[i*binFactor:(i+1)*binFactor])
	binz = np.zeros(nz//binFactor)
	for i in range(len(binz)):
		binz[i] = np.mean(z[i*binFactor:(i+1)*binFactor])
	plt.plot(binz, binCorrsNorm, 'bo')
	plt.xlabel('z/H');
	plt.axvline(x=-1.12, linestyle='--', color='k'); plt.axvline(x=1.12, linestyle='--', color='k');
	plt.title("Am=" + str(amDict[idString]))
	plt.ylabel('degree of correlation');
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();

	corrsNorm = corrsMax / corrsRand
	binFactor = 30
	binCorrsNorm = np.zeros(nz//binFactor)
	for i in range(len(binCorrsNorm)):
		binCorrsNorm[i] = np.mean(corrsNorm[i*binFactor:(i+1)*binFactor])
	binz = np.zeros(nz//binFactor)
	for i in range(len(binz)):
		binz[i] = np.mean(z[i*binFactor:(i+1)*binFactor])
	plt.plot(binz, binCorrsNorm, 'bo')
	plt.xlabel('z/H');
	plt.axvline(x=-1.12, linestyle='--', color='k'); plt.axvline(x=1.12, linestyle='--', color='k');
	plt.title("Am=" + str(amDict[idString]))
	plt.ylabel('degree of correlation');
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();


	############################################################################
	# MODELING
	############################################################################


	# setup
	distFromTrans = np.zeros_like(z)
	distFromTrans[:nz//2] = np.absolute(z[:nz//2] + zTrans)
	distFromTrans[nz//2:] = np.absolute(z[nz//2:] - zTrans)
	# rough density profile
	rhoInit = np.zeros_like(z)
	for i in range(len(z)):
		rhoInit[i] = 1.0 * np.exp(-np.square(z[i]))  
	rhoSS    = np.zeros_like(z)
	rhoSSLog = np.zeros_like(z)
	for i in range(len(z)):
		if   np.absolute(z[i]) < 2.0: 
			rhoSS[i] = rhoInit[i] 
		elif 2.0 < z[i]: 
			rhoSSLog[i] = np.log10(rhoInit[nz//4-1]) - ( (np.log10(rhoInit[nz//4-1])-np.log10(1.e-3))/(2.0) ) * (z[i]-2.0)  
			rhoSS[i]    = np.power(10, rhoSSLog[i]) 
		elif z[i] < -2.0: 
			rhoSSLog[i] = np.log10(rhoInit[nz//4-1]) + ( (np.log10(rhoInit[nz//4-1])-np.log10(1.e-3))/(2.0) ) * (z[i]+2.0)  
			rhoSS[i]    = np.power(10, rhoSSLog[i]) 


	# TURBULENT DIFFUSION ##########################################################
	# min beta induced by MRI	
	betaMin = np.zeros_like(z)
	for i in range(len(z)):
		if np.absolute(z[i]) <= zTrans: betaMin[i] = 50.0*np.power(amDict[idString],-1.2)
		#if np.absolute(z[i]) <= zTrans: betaMin[i] = 3.0
		else                          : betaMin[i] = 3.0
	betaMinMidList[nRun]  = np.mean(betaMin[nz//2-2:nz//2+3])
	# define alpha based on betaMin
	alpha   = 0.5*np.power(betaMin, -1.0)
	# define eta based on alpha (trivial for code units)
	etaTurb = alpha
	# local turbulent transport time (eta/h)		
	vTurbLoc = etaTurb
	# sum to get turbulent lag time
	lagTurb  = np.zeros_like(z)
	for i in range(len(z)):
		if   -100.0  < z[i] < -zTrans: lagTurb[i] = -np.sum(dx/vTurbLoc[i:ziTrans1]) 	
		elif -zTrans < z[i] < 0.0    : lagTurb[i] = -np.sum(dx/vTurbLoc[ziTrans1:i]) 
		elif 0.0     < z[i] < zTrans : lagTurb[i] = -np.sum(dx/vTurbLoc[i:ziTrans2])
		elif zTrans  < z[i] < 100.0  : lagTurb[i] = -np.sum(dx/vTurbLoc[ziTrans2:i])


	# AMBIPOLAR DIFFUSION	##########################################################
	# am is pre-defined
	am = np.zeros_like(z)
	for i in range(len(z)):
		if np.absolute(z[i]) < zTrans: am[i] = amDict[idString]
		else                         : am[i] = 10^20
	# actual beta
	byRmsProfile    = np.sqrt(np.mean(np.square(data[20:]), axis=(0,1,2)))
	pressureProfile = 0.5*rhoSS 
	betaActual      = pressureProfile / np.square(byRmsProfile)
	betaMidList[nRun] = np.mean(betaActual[nz//2-2:nz//2+3])
	betaTransList[nRun]  = 0.5 *( np.mean(betaActual[ziTrans1-2:ziTrans1+3]) + np.mean(betaActual[ziTrans2-2:ziTrans2+3]) )
	betaCoronaList[nRun] = 0.5 *( np.mean(betaActual[1*30:2*30]) + np.mean(betaActual[6*30:7*30]) )
	betaProfileList[nRun] = betaActual
	plt.semilogy(z, betaActual,      label='beta');
	plt.semilogy(z, betaMin,         label='betaMin');
	plt.semilogy(z, pressureProfile, label='pressure');
	plt.semilogy(z, byRmsProfile,    label='by');  
	plt.semilogy(z, alpha,           label='alpha');  
	plt.xlabel('z/H'); plt.ylabel('beta'); plt.legend(loc=(1.02, 0.0))
	plt.axvline(x=-1.12, linestyle='--', color='k'); plt.axvline(x=1.12, linestyle='--', color='k');
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();


	etaAd  = np.power(betaActual*am, -1.0)
	vAdLoc = etaAd
	lagAd  = np.zeros_like(z)
	for i in range(len(z)):
		if   -100.0  < z[i] < -zTrans: lagAd[i] = -np.sum(dx/vAdLoc[i:ziTrans1]) 	
		elif -zTrans < z[i] < 0.0    : lagAd[i] = -np.sum(dx/vAdLoc[ziTrans1:i]) 
		elif 0.0     < z[i] < zTrans : lagAd[i] = -np.sum(dx/vAdLoc[i:ziTrans2])
		elif zTrans  < z[i] < 100.0  : lagAd[i] = -np.sum(dx/vAdLoc[ziTrans2:i])


	# Buoyancy	##########################################################
	g              = np.absolute(z)
	aBuoyancyLoc   = g * np.power(betaActual, -1) * 1.0
	#vBuoyancyLoc   = np.power(vaLoc,1.5)
	vParticle      = np.zeros_like(z)
	lagBuoyancy    = np.zeros_like(z)
	timeSpent      = np.zeros_like(z)
	#for i in range(len(z)):
		#if   -100.0  < z[i] < -zTrans: lagBuoyancy[i] = -np.sum(dx/vBuoyancyLoc[i:ziTrans1]) 	
		#elif -zTrans < z[i] < 0.0    : lagBuoyancy[i] = 10.0 
		#elif 0.0     < z[i] < zTrans : lagBuoyancy[i] = 10.0
		#elif zTrans  < z[i] < 100.0  : lagBuoyancy[i] = -np.sum(dx/vBuoyancyLoc[ziTrans2:i])
	vParticle[ziTrans2] = 0.1
	for i in range(ziTrans2+1, nz):
		timeSpent[i-1] = dx/vParticle[i-1]
		vParticle[i]   = vParticle[i-1] + aBuoyancyLoc[i-1]*timeSpent[i-1] 
		lagBuoyancy[i] = lagBuoyancy[i-1] - timeSpent[i-1] 
			

	plt.plot(z, lagTurb, label='turb'); 
	plt.plot(z, lagAd,   label='AD'  ); 
	plt.plot(z, lagBuoyancy, label='buoyancy'  ); 

	plt.plot(z, lags, 'bo')

	plt.xlabel('z/H');
	plt.title("Am=" + str(amDict[idString]))
	plt.legend(loc=(1.02,0.0));
	plt.ylabel('most correlated time shift (orbits)');
	plt.axvline(x=-1.12, linestyle='--', color='k'); plt.axvline(x=1.12, linestyle='--', color='k');
	plt.ylim(-15,1); 
	plt.axhline(y=0, linestyle='--', color='k') 
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();


	plt.plot(lags[:nz//2], corrsMax[:nz//2], 'bo')
	plt.plot(lags[nz//2:], corrsMax[nz//2:], 'ro')
	plt.xlabel('lag of maximum correlation');
	plt.title("Am=" + str(amDict[idString]))
	plt.ylabel('absolute correlation');
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();

	plt.plot(lags[:nz//2], corrsNorm[:nz//2], 'bo')
	plt.plot(lags[nz//2:], corrsNorm[nz//2:], 'ro')
	plt.xlabel('lag of maximum correlation');
	plt.title("Am=" + str(amDict[idString]))
	plt.ylabel('normalized correlation');
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();


	corrsNorm[0] = 1.0
	corrsNorm[-1] = 1.0
	plt.errorbar(z, lags, yerr=2.0/corrsNorm)
	plt.xlabel('z/H'); plt.ylabel('lag of max correlation (orbits)')
	plt.title("Am=" + str(amDict[idString]))
	plt.axvline(x=-1.12, linestyle='--', color='k'); plt.axvline(x=1.12, linestyle='--', color='k');
	plt.axhline(y=0, linestyle='--', color='k') 
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();

	plt.plot(z, lags)
	plt.fill_between(z, lags-2.0/corrsNorm, lags+2.0/corrsNorm)
	plt.xlabel('z/H'); plt.ylabel('lag of max correlation (orbits)')
	plt.title("Am=" + str(amDict[idString]))
	plt.axvline(x=-1.12, linestyle='--', color='k'); plt.axvline(x=1.12, linestyle='--', color='k');
	plt.axhline(y=0, linestyle='--', color='k') 
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();


	amList       = [0.01, 0.018, 0.032, 0.056, 0.1 , 0.18, 0.32, 0.56, 1.0, 1.8 , 3.2,  5.6, 10.0 ]
	modelNorm    = 0.08
	modelLagList = [np.power(amList[i]/modelNorm, 0.33) for i in range(len(amList))]
	modelLagList[10] = modelLagList[9]
	modelLagList[11] = modelLagList[9]
	modelLagList[12] = modelLagList[9]
	plt.loglog(amList, midPlaneLagList, 'bo', label='mid-plane')
	plt.loglog(amList, modelLagList, label='by-eye power law')
	plt.loglog(amList, coronaLagList, 'go', label='corona')
	plt.axhline(y=1, linestyle='--', color='k');
	plt.axvline(x=modelNorm, linestyle='--', color='k');
	plt.axvline(x=2, linestyle='--', color='k');
	plt.xlabel('Am')
	plt.ylabel('lag')
	plt.legend();
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();

	plt.loglog(amList, betaMidList, 'bo', label='beta mid')
	plt.loglog(amList, betaMinMidList, 'r-', label='beta mid min MRI')
	plt.loglog(amList, betaTransList, 'ko', label='beta transition')
	plt.loglog(amList, betaCoronaList, 'go', label='beta corona')
	plt.axhline(y=3, linestyle='--', color='k')
	#plt.axhline(y=np.mean(betaMidList[np.nonzero(betaMidList)]), linestyle='--', color='b')
	#plt.axhline(y=np.mean(betaTransList[np.nonzero(betaTransList)]), linestyle='--', color='k')
	#plt.axhline(y=np.mean(betaCoronaList[np.nonzero(betaCoronaList)]), linestyle='--', color='g')
	plt.ylabel(r'$\beta_y$')
	plt.xlabel('Am')
	plt.legend()
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();

	for i in range(2,nRun+1):
		plt.semilogy(z, betaProfileList[i], label=str(amList[i]))
	plt.xlabel('z/H'); plt.ylabel(r"$\beta_y$"); 	plt.legend(loc=(1.02, 0.0));
	plt.axvline(x=-zTrans, linestyle='--', color='k'); plt.axvline(x=zTrans, linestyle='--', color='k');
	plt.axhline(y=3, linestyle='--', color='k'); 
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();

	plt.loglog(amList, midCorrList, 'bo')
	plt.ylabel('absolute correlation mid-plane to transition')
	plt.savefig(pdf, format='pdf', bbox_inches='tight'); plt.clf();
	

	# 4102 4103 4104 4105 4106 4107 4108 4109 4110 4111 4112



	pdf.close()













