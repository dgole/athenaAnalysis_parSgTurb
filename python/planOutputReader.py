#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors
import sys
sys.path.append('../python')
import athenaTools as tools
import numpy.polynomial.polynomial as poly
import emcee
import time
################################################################################
# Data class ###################################################################
################################################################################
class DataPlan:
	def __init__(self, path, dt=0.1, nPar=128*128*128, G=0.1, nStart=200, nTot=500, tSg=20.0):
		print("initializing PLAN data structure from " + path)

		# housekeeping
		self.strBase  = 'peaks_at_'
		self.strEnd   = '.txt'
		self.nStart   = nStart
		self.nTot     = nTot
		self.dt       = dt
		self.path     = path
		self.nPar     = nPar
		self.G        = G
		self.m0_ceres = G * (720.0/0.05)
		self.timeList = [n*self.dt for n in range(self.nTot)]
		self.tMax     = self.dt * nTot - tSg
		self.mg       = 9.766e-6

		# read in each peak file and save as an array in a list
		self.peakArrayList = [None for n in range(self.nTot)]
		self.nClumpsList   = [0    for n in range(self.nTot)]
		for n in range(self.nStart, self.nTot):
			# read in propper peak file
			print('###################################################')
			print('looking for peak file for n=' + str(n))
			t  = n*self.dt
			n1 = int(np.floor(t))
			n2 = np.int(np.floor(n-(1.0/self.dt)*n1))
			if   len(str(n1)) == 3: str1=str(n1)
			elif len(str(n1)) == 2: str1='0'+str(n1)
			elif len(str(n1)) == 1: str1='00'+str(n1)
			try:
				str2 = str(n2)+'00'
				fileName = self.strBase+str1+'.'+str2+self.strEnd
				temp     = np.loadtxt(self.path+fileName)
			except:
				str2 = str(n2)+'01'
				fileName = self.strBase+str1+'.'+str2+self.strEnd
				temp     = np.loadtxt(self.path+fileName)
			print('read in file named ' + fileName)
			# assign data to array / list
			if len(temp.shape) == 2:
				self.peakArrayList[n] = temp
			elif len(temp.shape) == 1 and temp.shape[0] == 19:
				self.peakArrayList[n] = temp.reshape((1, 19))
			try:
				self.nClumpsList[n] = (self.peakArrayList[n]).shape[0]
				print(str(self.nClumpsList[n]) + ' clumps found')
				#self.peakArrayList[n][:,2] *= self.m0_ceres
				self.peakArrayList[n][:,2] /= self.mg
			except:
				print('no clumps found')

		self.nClumps = np.asarray(self.nClumpsList)
		self.time    = np.asarray(self.timeList) - tSg

		if np.sum(self.nClumps)==0: print('NO CLUMPS AT ANY TIME!')
		mTestPar = self.peakArrayList[-1][0,2]
		nTestPar = self.peakArrayList[-1][0,1]
		self.mParTot = self.nPar * (mTestPar/nTestPar)


def getMasses(do, n):
	if do.nClumps[n] != 0:
		masses = []
		for mass in do.peakArrayList[n][:,2]:
			masses.append(mass)
		masses.sort()
		return masses

def getDiffMassHist(do, n):
	masses = getMasses(do, n)
	if masses is not None:
		dndmpList  = []
		mpList     = []
		for i in range(1, len(masses)-1):
			mp = masses[i]
			mpList.append(mp)
			dndmpList.append(2.0/(masses[i+1] - masses[i-1]))
		return mpList, dndmpList

def getDiffMassHist2(masses):
	if masses is not None:
		dndmpList  = []
		mpList     = []
		for i in range(1, len(masses)-1):
			mp = masses[i]
			mpList.append(mp)
			dndmpList.append(2.0/(masses[i+1] - masses[i-1]))
		return mpList, dndmpList


def getCumMassHist(do, n):
	if do.nClumps[n] != 0:
		masses  = getMasses(do, n)
		masses  = np.asarray(masses)
		nMasses = masses.shape[0]
		ngtm   = np.zeros_like(masses)
		for n in range(nMasses):
			ngtm[n] = nMasses - n
		return masses, ngtm
	else:
		print('no partilce found in this output')

def getCumMassHist2(masses):
	masses.sort()
	masses  = np.asarray(masses)
	nMasses = masses.shape[0]
	ngtm   = np.zeros_like(masses)
	for n in range(nMasses):
		ngtm[n] = nMasses - n
	return masses, ngtm

################################################################################
# modeling
################################################################################

def get_p_mle(do, n):
	nplan   = do.nClumpsList[n]
	if nplan > 3:
		minMass = np.amin(do.peakArrayList[n][:,2])
		sum     = 0
		for mass in do.peakArrayList[n][:,2]:
			sum += np.log(mass / minMass)
		p   = 1 + nplan * np.power(sum, -1)
		err = (p-1)/np.sqrt(nplan)
		return p, err

def get_p_mle2(sortedMassArr):
	minMass = np.amin(sortedMassArr)
	nplan   = sortedMassArr.shape[0]
	sum     = 0
	for mass in sortedMassArr:
		sum += np.log(mass / minMass)
	p   = 1 + nplan * np.power(sum, -1)
	err = (p-1)/np.sqrt(nplan)
	return p, err

def get_p_fit(do, n):
	nplan   = do.nClumpsList[n]
	if nplan > 3:
		mp, dndmp = getDiffMassHist(do, n)
		mp = np.asarray(mp); dndmp = np.asarray(dndmp);
		x = np.log10(mp); y = np.log10(dndmp);
		coefs, stats = poly.polyfit(x, y, 1, full=True)
		this_p = -coefs[1]; this_c = np.power(10,coefs[0])
		return this_p, 0.0, this_c

def convert_to_x(masses):
	return np.log(masses/np.amin(masses))

def p_to_P_bpl(masses, p):
	maxMass = np.amax(masses)
	elongateMasses = np.logspace(np.log10(maxMass*1.01), np.log10(maxMass*1000), num=1000)
	longMasses = np.append(masses, elongateMasses)

	nmOG   = masses.shape[0]
	nmLong = longMasses.shape[0]
	dm = np.zeros(nmLong)
	for i in range(1, nmLong-1):
		dm[i] = (longMasses[i+1] - longMasses[i-1])/ 2.0
	dm[0]  = dm[1]
	dm[-1] = dm[-2]
	dm_times_p = dm*p
	P  = np.zeros_like(p)
	for i in range(P.shape[0]):
		P[i] = np.sum(dm_times_p[i:])
	return P


	'''
	dlnm  = np.zeros_like(masses)
	for n in range(masses.shape[0]-1):
		dlnm[n] = np.log(masses[n+1]/masses[n])
	dlnm[-1]=dlnm[-2]
	dm    = np.zeros_like(masses)
	for n in range(masses.shape[0]-1):
		dm[n] = masses[n+1]-masses[n]
	dm[-1]=dm[-2]
	dndmp = pdf * dlnm/dm * nRealMasses
	return dndmp
	'''



def gridSearch1d(funcName, mp, p1min, p1max, p1step):
	#print("performing grid search...")
	maxLike = -1.e10
	p1 = np.arange(p1min, p1max, p1step)
	grid = np.zeros(p1.shape[0])
	for n1 in range(p1.shape[0]):
		params      = (p1[n1])
		lnLike      = funcName(mp, params)
		grid[n1] = lnLike
		if lnLike > maxLike:
			maxLike     = lnLike
			paramsOfMax = params
			#print("#######################################")
			#print("new max found: " + str(maxLike))
			#print("params:        " + str(params))
	return paramsOfMax, maxLike

def gridSearch2d(funcName, mp, p1min, p1max, p1step, p2min, p2max, p2step):
	#print("performing grid search...")
	maxLike = -1.e10
	p1 = np.arange(p1min, p1max, p1step)
	p2 = np.arange(p2min, p2max, p2step)
	grid = np.zeros((p1.shape[0], p2.shape[0]))
	for n1 in range(p1.shape[0]):
		for n2 in range(p2.shape[0]):
			params      = (p1[n1], p2[n2])
			lnLike      = funcName(mp, params)
			grid[n1,n2] = lnLike
			if lnLike > maxLike:
				maxLike     = lnLike
				paramsOfMax = params
				#print("#######################################")
				#print("new max found: " + str(maxLike))
				#print("params:        " + str(params))
	return paramsOfMax, maxLike

def gridSearch3d(funcName, mp, p1min, p1max, p1step, p2min, p2max, p2step, p3min, p3max, p3step):
	#print("performing grid search...")
	maxLike = -1.e10
	p1 = np.arange(p1min, p1max, p1step)
	p2 = np.arange(p2min, p2max, p2step)
	p3 = np.arange(p3min, p3max, p3step)
	grid = np.zeros((p1.shape[0], p2.shape[0], p3.shape[0]))
	for n1 in range(p1.shape[0]):
		for n2 in range(p2.shape[0]):
			for n3 in range(p3.shape[0]):
				params         = (p1[n1], p2[n2], p3[n3])
				lnLike         = funcName(mp, params)
				grid[n1,n2,n3] = lnLike
				if lnLike > maxLike:
					maxLike     = lnLike
					paramsOfMax = params
					#print("#######################################")
					#print("new max found: " + str(maxLike))
					#print("params:        " + str(params))
	return paramsOfMax, maxLike

def bootstrap(mp, funcName, nParams, nb=100, sampleFactor=1.0, pathSave=None):
	print("Computing " + str(funcName.__name__)
		+ " fit for " + str(nb) + " subsamples...")
	paramsAll  =  np.zeros([nb, nParams])
	maxLikeAll = np.zeros(nb)
	t0 = time.time()
	for n in range(nb):
		if n%20==0 and n>0:
			str1 = "n = " + str(n) + " of " + str(nb)
			str2 = "average time per sample = " + str(np.round((time.time()-t0)/n,3))
			print(str1 + ", " + str2)
		sample          = np.random.choice(mp,size=int(sampleFactor*mp.shape[0]),replace=True)
		params, maxLike = funcName(sample)
		paramsAll[n]    = params
	means           = np.median(paramsAll, axis=0)
	params, maxLike = funcName(mp)
	errsPlus        = np.percentile(paramsAll, 84, axis=0) - params
	errsMinus       = params - np.percentile(paramsAll, 14, axis=0)
	#print(params)
	#print(np.percentile(paramsAll, 86, axis=0))
	#print(np.percentile(paramsAll, 14, axis=0))
	#print(errsPlus)
	#print(errsMinus)
	if pathSave is not None:
		for i in range(paramsAll.shape[1]):
			plt.figure(0)
			nBins = int(max(10.0, nb/50.0))
			plt.hist(paramsAll[:,i], color='gray', bins=nBins)
			plt.xlim(np.amin(paramsAll[:,i])-0.1, np.amax(paramsAll[:,i])+0.1)
			try:
				plt.axvline(params[i], color='k', linestyle='--')
				plt.axvline(params[i]+errsPlus[i], color='r', linestyle='--')
				plt.axvline(params[i]-errsMinus[i], color='b', linestyle='--')
				plt.title(str(np.round(params[i],4)) + " + " + str(np.round(errsPlus[i],4)) + " - " + str(np.round(errsMinus[i],4)))
			except:
				plt.axvline(params, color='k', linestyle='--')
				plt.axvline(params+errsPlus, color='r', linestyle='--')
				plt.axvline(params-errsMinus, color='b', linestyle='--')
				plt.title(str(np.round(params,4)) + " + " + str(np.round(errsPlus[i],4)) + " - " + str(np.round(errsMinus[i],4)))
			tools.saveAndClear(pathSave + "z_bs_" + str(funcName.__name__) + "_" + str(i) +".png", figNum=0)
	return params, errsPlus, errsMinus, maxLike




# SPL FITTING
def p_spl(masses, params):
	alpha  = params
	xArr   = convert_to_x(masses)
	return alpha * np.exp(-alpha*xArr)
def P_spl(masses, params):
	alpha   = params
	minMass = np.amin(masses)
	nMasses = masses.shape[0]
	ngtm    = np.power(masses/minMass, -alpha)
	return ngtm
def lnlike_spl(masses, params):
	pArr = p_spl(masses, params)
	sum  = np.sum(np.log(pArr))
	return sum
def fit_spl(mp):
	a = 5.0
	minMass = np.amin(mp); maxMass = np.amax(mp)
	p1ll = -5.0; p1ul = 5.0; p1step = 0.5;
	p1, grid1 = gridSearch1d(lnlike_spl, mp,
							 p1ll, p1ul, p1step)
	for i in range(5):
		p1, grid2 = gridSearch1d(lnlike_spl, mp,
		 						 max(p1-2*p1step,p1ll), min(p1+2*p1step,p1ul), p1step/a)
		p1step /= a;
	return p1, grid1




# STPL FITTING
def p_stpl(masses, params):
	alpha, x_exp = params
	xArr   = convert_to_x(masses)
	top    = alpha + np.exp(-x_exp) * np.exp(xArr)
	bottom = np.exp( alpha*xArr + np.exp(-x_exp)*(np.exp(xArr)-1.0) )
	return top/bottom
def P_stpl(masses, params):
	alpha, x_exp = params
	minMass = np.amin(masses)
	mExp    = np.exp(x_exp)*minMass
	nMasses = masses.shape[0]
	ngtm    = np.power(masses/minMass, -alpha) * np.exp(-(masses-minMass)/mExp)
	return ngtm
def lnlike_stpl(masses, params):
	pArr = p_stpl(masses, params)
	sum  = np.sum(np.log(pArr))
	return sum
def fit_stpl(mp):
	a = 5.0
	minMass = np.amin(mp); maxMass = np.amax(mp)
	p1ll = -5.0; p1ul = 5.0; p1step = 0.5;
	p2ll = 0.01;  p2ul = np.log(maxMass/minMass); p2step = 0.5;
	(p1, p2), grid1 = gridSearch2d(lnlike_stpl, mp,
								   p1ll, p1ul, p1step,
								   p2ll, p2ul, p2step)
	for i in range(5):
		(p1, p2), grid2 = gridSearch2d(lnlike_stpl, mp,
		 							   max(p1-2*p1step,p1ll), min(p1+2*p1step,p1ul), p1step/a,
		 							   max(p2-2*p2step,p2ll), min(p2+2*p2step,p2ul), p2step/a)
		p1step /= a; p2step /= a;
	return (p1, p2), grid1




# VTPL FITTING
def p_vtpl(masses, params):
	alpha, beta, x_exp = params
	xArr   = convert_to_x(masses)
	top    = alpha + beta*np.exp(beta*(x_exp +xArr))
	bottom = np.exp( alpha*xArr + np.exp(beta*x_exp)*(np.exp(beta*xArr)-1.0) )
	return top/bottom
def P_vtpl(masses, params):
	alpha, beta, x_exp = params
	minMass = np.amin(masses)
	mExp    = np.exp(x_exp)*minMass
	nMasses = masses.shape[0]
	part1   = np.power(masses/minMass, -alpha)
	part2   = np.exp(-(np.power(masses,beta)-np.power(minMass,beta))/np.power(mExp,beta))
	ngtm    = part1*part2
	return ngtm
def lnlike_vtpl(masses, params):
	pArr = p_vtpl(masses, params)
	sum  = np.sum(np.log(pArr))
	return sum
def fit_vtpl(mp):
	a = 5.0
	minMass = np.amin(mp); maxMass = np.amax(mp)
	p1ll = -5.0; p1ul = 5.0; p1step = 0.5;
	p2ll = -5.0; p2ul = 5.0; p2step = 0.5;
	p3ll = 0.01;  p3ul = np.log(maxMass/minMass); p3step = 0.5;
	(p1, p2, p3), grid1 = gridSearch3d(lnlike_vtpl, mp,
									  p1ll, p1ul, p1step,
									  p2ll, p2ul, p2step,
									  p3ll, p3ul, p3step)
	for i in range(5):
		#print("iteration: " + str(i))
		(p1, p2, p3), grid2 = gridSearch3d(lnlike_vtpl, mp,
	 										max(p1-2*p1step,p1ll), min(p1+2*p1step,p1ul), p1step/a,
	 										max(p2-2*p2step,p2ll), min(p2+2*p2step,p2ul), p2step/a,
	 										max(p3-2*p3step,p3ll), min(p3+2*p3step,p3ul), p3step/a)
		#print(p1, p2, p3)
		#print(p1step, p2step, p3step)
		p1step /= a; p2step /= a; p3step /= a;
	return (p1, p2, p3), grid1




# BCPL FITTING
def p_bcpl(masses, params):
	a1, a2, xb = params
	xArr   = convert_to_x(masses)
	ltxbr  = np.where(xArr <= xb, 1.0, 0.0)
	gtxbr  = np.where(xArr >  xb, 1.0, 0.0)
	ltxbr *= a1 * np.exp(-a1 * xArr)
	gtxbr *= a2 * np.exp((a2-a1)*xb - a2*xArr)
	return gtxbr + ltxbr
def P_bcpl(masses, params):
	a1, a2, xb = params
	minMass = np.amin(masses)
	mb      = np.exp(xb)*minMass
	nMasses = masses.shape[0]
	ltxbr  = np.where(masses <= mb, 1.0, 0.0)
	gtxbr  = np.where(masses >  mb, 1.0, 0.0)
	ltxbr *= np.power(masses/minMass, -a1)
	gtxbr *= np.power(masses, -a2)/(np.power(mb,a1-a2)*np.power(minMass,-a1))
	ngtm    = (ltxbr+gtxbr)
	return ngtm
def lnlike_bcpl(masses, params):
	pArr = p_bcpl(masses, params)
	sum  = np.sum(np.log(pArr))
	return sum
def fit_bcpl(mp):
	a = 5.0
	minMass = np.amin(mp); maxMass = np.amax(mp)
	p1ll = -5.0; p1ul = 5.0; p1step = 0.5;
	p2ll = 0.0;  p2ul = 5.0; p2step = 0.5;
	p3ll = 0.01;  p3ul = np.log(maxMass/minMass); p3step = 0.5;
	(p1, p2, p3), grid1 = gridSearch3d(lnlike_bcpl, mp,
									  p1ll, p1ul, p1step,
									  p2ll, p2ul, p2step,
									  p3ll, p3ul, p3step)
	for i in range(5):
		(p1, p2, p3), grid2 = gridSearch3d(lnlike_bcpl, mp,
	 										max(p1-2*p1step,p1ll), min(p1+2*p1step,p1ul), p1step/a,
	 										max(p2-2*p2step,p2ll), min(p2+2*p2step,p2ul), p2step/a,
	 										max(p3-2*p3step,p3ll), min(p3+2*p3step,p3ul), p3step/a)
		p1step /= a; p2step /= a; p3step /= a;
	return (p1, p2, p3), grid1




# TPL FITTING
def p_tpl(masses, params):
	a, xt  = params
	xArr   = convert_to_x(masses)
	top    = a * np.exp(-a*xArr)
	bottom = 1.0 - np.exp(-a*xt)
	return top/bottom
def P_tpl(masses, params):
	a, x_tr = params
	minMass = np.amin(masses)
	mTr     = np.exp(x_tr)*minMass
	nMasses = masses.shape[0]
	t1      = np.power(masses/minMass, -a)
	t2      = np.power(mTr/minMass, -a)
	t3      = np.power(mTr/minMass, -a)
	ngtm    = (t1-t2)/(1-t3)
	return ngtm
def lnlike_tpl(masses, params):
	pArr = p_tpl(masses, params)
	sum  = np.sum(np.log(pArr))
	return sum
def fit_tpl(mp):
	a=5.0
	minMass = np.amin(mp); maxMass = np.amax(mp)
	#print(maxMass)
	p1ll = -5.0; p1ul = 5.0; p1step = 0.1;
	p2ll = np.log(maxMass/minMass);  p2ul = p2ll*10.0; p2step = (p2ul-p1ul)/50.0;
	(p1, p2), grid1 = gridSearch2d(lnlike_tpl, mp,
								   p1ll, p1ul, p1step,
								   p2ll, p2ul, p2step)
	for i in range(5):
		(p1, p2), grid2 = gridSearch2d(lnlike_tpl, mp,
		 							   max(p1-2*p1step,p1ll), min(p1+2*p1step,p1ul), p1step/a,
		 							   max(p2-2*p2step,p2ll), min(p2+2*p2step,p2ul), p2step/a)
		p1step /= a; p2step /= a;
	return (p1, p2), grid1




# BPL FITTING
def p_bpl(masses, params):
	a1, a2, xb = params
	c0     = np.power((1./a1)+((1./a2)-(1./a1))*np.power(np.exp(xb),-a1),-1)
	xArr   = convert_to_x(masses)
	ltxbr  = np.where(xArr <= xb, 1.0, 0.0)
	gtxbr  = np.where(xArr >  xb, 1.0, 0.0)
	ltxbr *= c0 * np.exp(-a1 * xArr)
	gtxbr *= c0 * np.exp((a2-a1)*xb - a2*xArr)
	return gtxbr + ltxbr
def P_bpl(masses, params):
	'''
	maxMass        = np.amax(masses)
	elongateMasses = np.logspace(np.log10(maxMass*1.01), np.log10(maxMass*1000), num=1000)
	longMasses     = np.append(masses, elongateMasses)
	p              = p_bpl(longMasses, params)
	nmOG           = masses.shape[0]
	nmLong         = longMasses.shape[0]
	dm             = np.zeros(nmLong)
	for i in range(1, nmLong-1):
		dm[i] = (longMasses[i+1] - longMasses[i-1])/ 2.0
	dm[0]  = dm[1]
	dm[-1] = dm[-2]
	dm_times_p = dm*p
	P  = np.zeros_like(p)
	for i in range(P.shape[0]):
		P[i] = np.sum(dm_times_p[i:])
	#plt.figure(1); plt.loglog(longMasses, p); plt.show(1); plt.clf(); plt.figure(0);
	#plt.figure(1); plt.loglog(longMasses, P); plt.show(1); plt.clf(); plt.figure(0);
	return P[:nmOG]
	'''

	maxMass        = np.amax(masses)
	elongateMasses = np.logspace(np.log10(maxMass*1.01), np.log10(maxMass*1000), num=100)
	longMasses     = np.append(masses, elongateMasses)
	p        = p_bpl(longMasses, params)
	xArr = convert_to_x(longMasses)
	dx   = np.zeros_like(longMasses)
	for i in range(len(dx)-1):
		dx[i] = xArr[i+1] - xArr[i]
	dx[-1] = dx[-2]
	print(dx)
	ngtm    = np.zeros_like(dx)
	for i in range(0, len(dx)):
		ngtm[i] = np.sum(p[i:-1]*dx[i:-1])
	return ngtm[:masses.shape[0]]


	'''
	p    = p_bpl(masses, params)
	xArr = convert_to_x(masses)
	dx   = np.zeros_like(masses)
	for i in range(len(dx)-1):
		dx[i] = xArr[i+1] - xArr[i]
	dx[-1] = dx[-2]
	ngtm    = np.zeros_like(dx)
	for i in range(0, len(dx)):
		ngtm[i] = np.sum(p[i:-1]*dx[i:-1])
	return ngtm
	'''

def lnlike_bpl(masses, params):
	pArr = p_bpl(masses, params)
	sum  = np.sum(np.log(pArr))
	return sum
def fit_bpl(mp):
	a = 5.0
	minMass = np.amin(mp); maxMass = np.amax(mp)
	p1ll = -5.0; p1ul = 5.0; p1step = 0.5;
	p2ll = 0.1;  p2ul = 5.0; p2step = 0.5;
	p3ll = 0.1;  p3ul = np.log(maxMass/minMass); p3step = 0.5;
	(p1, p2, p3), grid1 = gridSearch3d(lnlike_bpl, mp,
									  p1ll, p1ul, p1step,
									  p2ll, p2ul, p2step,
									  p3ll, p3ul, p3step)
	for i in range(5):
		(p1, p2, p3), grid2 = gridSearch3d(lnlike_bpl, mp,
	 										max(p1-2*p1step,p1ll), min(p1+2*p1step,p1ul), p1step/a,
	 										max(p2-2*p2step,p2ll), min(p2+2*p2step,p2ul), p2step/a,
	 										max(p3-2*p3step,p3ll), min(p3+2*p3step,p3ul), p3step/a)
		#print("iteration " + str(i))
		#print(p1, p2, p3)
		#print(max(p1-2*p1step,p1ll), min(p1+2*p1step,p1ul), p1step/a)
		#print(max(p2-2*p2step,p2ll), min(p2+2*p2step,p2ul), p2step/a)
		#print(max(p3-2*p3step,p3ll), min(p3+2*p3step,p3ul), p3step/a)
		p1step /= a; p2step /= a; p3step /= a;
	return (p1, p2, p3), grid1



def BIC(K, N, lnLike):
	return 2.0*K*np.log(N)-2.0*lnLike
def AIC(K, N, lnLike):
	return 2.0*K-2.0*lnLike

def reportParams(name, paramNames, means, errsPlus, errsMinus, maxLike, dbic, daic, mp1_min):
	nParams = len(paramNames)
	print("########################################################")
	print(name)
	if nParams == 1:
		print(paramNames[0].ljust(8) + " = "
			+ str(np.round(means,4)) + " p "
			+ str(np.round(errsPlus[0],4)) + " m "
			+ str(np.round(errsMinus[0],4)))
	else:
		for n in range(nParams):
			print(paramNames[n].ljust(8) + " = "
				+ str(np.round(means[n],4)) + " p "
				+ str(np.round(errsPlus[n],4)) + " m "
				+ str(np.round(errsMinus[n],4)))
			if paramNames[n][0]=='x':
				mExp         = mp1_min*np.exp(means[n])
				mExpErrPlus  = mp1_min*np.exp(means[n]+errsPlus[n]) - mExp
				mExpErrMinus = mExp - mp1_min*np.exp(means[n]-errsMinus[n])
				print(("M_"+paramNames[n][1:]).ljust(8) + " = "
					+ str(np.round(mExp,6)) + " p "
					+ str(np.round(mExpErrPlus,6)) + " m "
					+ str(np.round(mExpErrMinus,6)))
	print("ln(like) = " + str(np.round(maxLike,3)))
	print("DBIC     = " + str(np.round(dbic,3)))
	print("DAIC     = " + str(np.round(daic,3)))

####################################################################
# plotting functions ###############################################
####################################################################
def plotCumMassHist(do, nStart=None, nEnd=None, spacing=0.05, figNum=0, legendLabel=None, colorOption='ko'):
	print('plotting time-averaged cumulative mass hist for ' + do.path)
	plt.figure(figNum)
	if nStart is None: nStart = do.nt // 2
	if nEnd   is None: nEnd   = do.nt
	count=0
	for n in range(nStart, nEnd):
		bins, hist = getCumMassHist(do, n, spacing=spacing)
		if not isinstance(bins, int):
			if count==0:
				masterHist  = hist
				plotBins    = bins
			else:
				masterHist += hist
			count+=1
	masterHist /= count
	plt.loglog(plotBins, masterHist, colorOption, label=legendLabel)
	plt.xlabel(r'$M_p$')
	plt.ylabel(r'$N(>M_p)$')

def scatterPlotXZ(do, n, figNum=0):
	print('plotting XZ scatter plot for ' + do.path + ', n=' +  str(n))
	plt.figure(figNum)
	masses = do.peakArrayList[n][:, 2]
	xs     = do.peakArrayList[n][:, 4]
	ys     = do.peakArrayList[n][:, 5]
	zs     = do.peakArrayList[n][:, 6]
	sizes  = [np.power(3.e3*mass, 1./2.) for mass in masses]
	############################################################################
	plt.scatter(xs, zs, s=sizes)
	plt.xlabel(r'x')
	plt.ylabel(r'z')
	plt.ylim(-0.1, 0.1)
	plt.xlim(-0.1, 0.1)
	plt.title(r'$t=$' + str(np.round(do.timeList[n],1)))

def scatterPlotXY(do, n, figNum=0):
	print('plotting XY scatter plot for ' + do.path + ', n=' +  str(n))
	plt.figure(figNum)
	masses = do.peakArrayList[n][:, 2]
	xs     = do.peakArrayList[n][:, 4]
	ys     = do.peakArrayList[n][:, 5]
	zs     = do.peakArrayList[n][:, 6]
	sizes  = [np.power(3.e3*mass, 1./2.) for mass in masses]
	############################################################################
	plt.scatter(xs, ys, s=sizes)
	plt.xlabel(r'x')
	plt.ylabel(r'y')
	plt.ylim(-0.1, 0.1)
	plt.xlim(-0.1, 0.1)
	plt.title(r'$t=$' + str(np.round(do.timeList[n],1)))

def scatterPlotXYZ(do, n, figNum=0):
	print('plotting XYZ scatter plot for ' + do.path + ', n=' +  str(n))
	plt.figure(figNum)
	masses = do.peakArrayList[n][:, 2]
	xs     = do.peakArrayList[n][:, 4]
	ys     = do.peakArrayList[n][:, 5]
	zs     = do.peakArrayList[n][:, 6]
	sizes  = [np.power(3.e3*mass, 1./2.) for mass in masses]
	fig = plt.figure(figNum)
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(ys, xs, zs, s=sizes)
	ax.set_xlabel(r'$\phi$')
	ax.set_ylabel(r'$r$')
	ax.set_zlabel(r'$z$')
	ax.set_xlim(-0.1, 0.1)
	ax.set_ylim(-0.1, 0.1)
	ax.set_zlim(-0.1, 0.1)
	plt.title(r'$t=$' + str(np.round(do.timeList[n],1)))




def get_p_avg(do, nStart=None, nEnd=None):
	if nStart is None: nStart = do.nt-100
	if nEnd   is None: nEnd   = do.nt
	pList   = []
	errList = []
	for n in range(nStart, nEnd):
		p, err = get_p(do, n)
		pList.append(p)
		errList.append(err)
	pArr   = np.asarray(pList)
	errArr = np.asarray(errList)
	pAvg   = np.round(np.mean(pArr),   decimals=4)
	errAvg = np.round(np.mean(errArr), decimals=4)
	return pAvg, errAvg

def pValuePlot(do, nStart=None, nEnd=None, figNum=0):
	plt.figure(figNum)
	pList   = []
	errList = []
	for n in range(do.nFirstClump, do.nt):
		p, err = get_p(do, n)
		pList.append(p)
		errList.append(err)
	pArr   = np.asarray(pList)
	errArr = np.asarray(errList)
	plt.plot(do.timeList[do.nFirstClump:], pArr)
	plt.plot(do.timeList[do.nFirstClump:], pArr-errArr)
	plt.plot(do.timeList[do.nFirstClump:], pArr+errArr)
	plt.ylim(1.0, 3.0)
	pAvg, errAvg = get_p_avg(do, nStart=nStart, nEnd=nEnd)
	plt.ylabel(r'$p$')
	plt.xlabel(r'$t \Omega$')
	plt.title(r'$p_{end}=$' + str(pAvg) + r'$\pm$' + str(errAvg))



































#
