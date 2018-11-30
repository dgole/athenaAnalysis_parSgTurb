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

#m.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

####################################################################
# Data class ###############################################
####################################################################
class Data3d:
	def __init__(self, path, baseName="Par_Strat3d", dt=1.0):
		print("initializing 3d data structure from " + path)
		self.path    = path
		fileList     = os.listdir(self.path)
		nFiles       = len(fileList)
		self.nt      = nFiles
		self.names   = [self.path + baseName + "." + tools.getTimeStepString(n) + ".npy"
		 				for n in range(nFiles)]
		dataTemp     = np.load(self.names[0])
		self.cols    = {
						'rho'    : 3,
						'vx'     : 4,
						'vy'     : 5,
						'vz'     : 6,
						'dpar'   : 7,
						'vxpar'  : 8,
						'vypar'  : 9,
						'vypar'  : 10
						}
		self.funcs   = {
						'v'         : v,
						'KE'        : ke,
						'drho'      : drho,
						'dvx'       : dvx,
						'dvy'       : dvy,
						'dvz'       : dvz,
						'dv'        : dv,
						'drhoNorm'  : drhoNorm,
						'dvxNorm'   : dvxNorm,
						'dvyNorm'   : dvyNorm,
						'dvzNorm'   : dvzNorm,
						'dvNorm'    : dvNorm,
						'KEx'       : kex,
						'KEy'       : key,
						'KEz'       : kez,
						'rootKEx'   : rootKEx,
						'rootKEy'   : rootKEy,
						'rootKEz'   : rootKEz,
						'dKE'       : dKE,
						'rootRhoVx' : rootRhoVx,
						'rootRhoVy' : rootRhoVy,
						'rootRhoVz' : rootRhoVz,
						'rootRhoDvx': rootRhoVx,
						'rootRhoDvy': rootRhoDvy,
						'rootRhoDvz': rootRhoDvz
									         	}
		self.header  = {
						'v'         : r'$v$',
						'KE'        : r'$\frac{1}{2}\rho v^2$',
						'rho'       : r'$\rho$',
						'vx'        : r'$v_x$',
						'vy'        : r'$v_y$',
						'vz'        : r'$v_z$',
						'dpar'      : r'$dpar$',
						'vxpar'     : r'$v_x_{par}$',
						'vypar'     : r'$v_y_{par}$',
						'vypar'     : r'$v_z_{par}$',
						'drho'      : r'$\delta \rho$',
						'dvx'       : r'$\delta v_x$',
						'dvy'       : r'$\delta v_y$',
						'dvz'       : r'$\delta v_z$',
						'dv'        : r'$\delta v$',
						'drhoNorm'  : r'$\delta \rho / \rho$',
						'dvxNorm'   : r'$\delta v_x / v_x$',
						'dvyNorm'   : r'$\delta v_y / v_y$',
						'dvzNorm'   : r'$\delta v_z / v_z$',
						'dvNorm'    : r'$\delta v / v$',
						'KEx'       : r'$\frac{1}{2}\rho v_x^2$',
						'KEy'       : r'$\frac{1}{2}\rho v_y^2$',
						'KEz'       : r'$\frac{1}{2}\rho v_z^2$',
						'rootKEx'   : r'$\sqrt{\frac{1}{2}\rho v_x^2}$',
						'rootKEy'   : r'$\sqrt{\frac{1}{2}\rho v_y^2}$',
						'rootKEz'   : r'$\sqrt{\frac{1}{2}\rho v_z^2}$',
						'dKE'       : r'$\frac{1}{2}\rho \delta v^2$',
						'rootRhoVx' : r'$\sqrt{0.5\rho}v_x$',
						'rootRhoVy' : r'$\sqrt{0.5\rho}v_y$',
						'rootRhoVz' : r'$\sqrt{0.5\rho}v_z$',
						'rootRhoDvx' : r'$\sqrt{0.5\rho}\deltav_x$',
						'rootRhoDvy' : r'$\sqrt{0.5\rho}\deltav_y$',
						'rootRhoDvz' : r'$\sqrt{0.5\rho}\deltav_z$'
						}
		self.x      = dataTemp[:,0,0,0]
		self.y      = dataTemp[0,:,0,1]
		self.z      = dataTemp[0,0,:,2]
		self.nx     = self.x.shape[0]
		self.ny     = self.y.shape[0]
		self.nz     = self.z.shape[0]
		self.dx     = self.x[2] - self.x[1]
		#self.zmax   = np.round(-self.z[0],1)
		self.dt     = dt
		self.tmax   = dt*self.nt
		self.t      = np.arange(0, self.tmax-dt/2.0, dt)
		print("data structure initialized pointing to data of shape " + str(dataTemp.shape)
		      + ' x ' + str(self.nt) + ' time steps')
		del dataTemp
	def get3d(self, key, n, *args, **kwargs):
		if   key in self.cols:
			returnData = np.load(self.names[n])[:,:,:,self.cols[key]]
		elif key in self.funcs:
			returnData = self.funcs[key](self, n)
		else:
			print('unknown key given')
			returnData = None
		return returnData
	def getxindex(self, x):
		return (np.abs(self.x-x)).argmin()
	def getzindex(self, y):
		return (np.abs(self.y-y)).argmin()
	def getzindex(self, z):
		return (np.abs(self.z-z)).argmin()
	def gettindex(self, t):
		return (np.abs(self.t-t)).argmin()

#########################################################################
# Plots
#########################################################################

def profile(do, key, figNum=0, tStart=None, tEnd=None, legendLabel=None, absAvg=1, absPlot=1):
	print(do.path + ": making profile plot for key " + key)
	sys.stdout.flush()
	plt.figure(figNum)
	if tStart == None: tStart = do.t[-1]/2.0
	if tEnd   == None: tEnd   = do.t[-1]
	nStart   = do.gettindex(tStart)
	nEnd     = do.gettindex(tEnd)
	plotData = np.zeros(do.nz)
	nCount   = 0
	for n in range(nStart, nEnd):
		if absAvg==1: data3d = np.absolute(do.get3d(key, n))
		else:            data3d =             do.get3d(key, n)
		data1d    = np.mean(data3d, axis=(0,1))
		plotData += data1d
		nCount   += 1
	plotData /= nCount
	title = do.header[key]
	plt.ylim(0.5*np.amin(np.absolute(plotData)),2.0*np.amax(np.absolute(plotData)))
	if absPlot==1: plt.semilogy(do.z, np.absolute(plotData), label=legendLabel)
	else:          plt.plot    (do.z,             plotData , label=legendLabel)
	plt.ylabel(do.header[key]);
	plt.xlabel(r"$z/H$");
	plt.tight_layout()



#########################################################################
# Helper
#########################################################################

def getPert(a):
	aAvg1d = np.mean(a, axis=(0,1))
	aAvg3d = np.zeros_like(a)
	for zi in range(a.shape[2]): aAvg3d[:,:,zi] = aAvg1d[zi]
	da     =  a - aAvg3d
	return da

def getPertNorm(a):
	aAvg1d = np.mean(a, axis=(0,1))
	aAvg3d = np.zeros_like(a)
	for zi in range(a.shape[2]): aAvg3d[:,:,zi] = aAvg1d[zi]
	daNorm = (a - aAvg3d) / aAvg3d
	return daNorm


#########################################################################
# Callable with do.get3d()
#########################################################################

def drho(do, n):
	print('calculating drho for n=' + str(n))
	sys.stdout.flush()
	da = getPert(do.get3d('rho', n))
	return da

def drhoAbs(do, n):
	print('calculating drhoAbs for n=' + str(n))
	sys.stdout.flush()
	da = np.absolute(getPert(do.get3d('rho', n)))
	return da

def drhoNorm(do, n):
	print('calculating drhoNorm for n=' + str(n))
	sys.stdout.flush()
	daNorm = getPertNorm(do.get3d('rho', n))
	return daNorm

def dvx(do, n):
	print('calculating dvx for n=' + str(n))
	sys.stdout.flush()
	da = getPert(do.get3d('vx', n))
	return da

def dvxNorm(do, n):
	print('calculating dvxNorm for n=' + str(n))
	sys.stdout.flush()
	daNorm = getPertNorm(do.get3d('vx', n))
	return daNorm

def dvy(do, n):
	print('calculating dvy for n=' + str(n))
	sys.stdout.flush()
	da = getPert(do.get3d('vy', n))
	return da

def dvyNorm(do, n):
	print('calculating dvyNorm for n=' + str(n))
	sys.stdout.flush()
	daNorm = getPertNorm(do.get3d('vy', n))
	return daNorm

def dvz(do, n):
	print('calculating dvz for n=' + str(n))
	sys.stdout.flush()
	da = getPert(do.get3d('vz', n))
	return da

def dvzNorm(do, n):
	print('calculating dvzNorm for n=' + str(n))
	sys.stdout.flush()
	daNorm = getPertNorm(do.get3d('vz', n))
	return daNorm

def dv(do, n):
	print('calculating dv for n=' + str(n))
	sys.stdout.flush()
	dvx = do.get3d('dvx', n)
	dvy = do.get3d('dvy', n)
	dvz = do.get3d('dvz', n)
	dv  = np.sqrt(np.square(dvx) + np.square(dvy) + np.square(dvz))
	return dv

def dvNorm(do, n):
	print('calculating dvNorm for n=' + str(n))
	sys.stdout.flush()
	dvx = do.get3d('dvx', n)
	dvy = do.get3d('dvy', n)
	dvz = do.get3d('dvz', n)
	dv  = np.sqrt(np.square(dvx) + np.square(dvy) + np.square(dvz))
	v   = np.sqrt(np.square(do.get3d('vx',n)) + np.square(do.get3d('vy',n)) + np.square(do.get3d('vz',n)))
	dvNorm = dv/v
	return dvNorm

def v(do, n):
	print('calculating v for n=' + str(n))
	sys.stdout.flush()
	vx = do.get3d('vx', n)
	vy = do.get3d('vy', n)
	vz = do.get3d('vz', n)
	v   = np.sqrt(np.square(vx) + np.square(vy) + np.square(vz))
	return v

def ke(do, n):
	print('calculating KE for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('v', n)
	rho = do.get3d('rho', n)
	ke  = 0.5*rho*np.square(v)
	return ke

def kex(do, n):
	print('calculating KEx for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('vx', n)
	rho = do.get3d('rho', n)
	return 0.5*rho*np.square(v)

def key(do, n):
	print('calculating KEy for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('vx', n)
	rho = do.get3d('rho', n)
	return 0.5*rho*np.square(v)

def kez(do, n):
	print('calculating KEz for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('vx', n)
	rho = do.get3d('rho', n)
	return 0.5*rho*np.square(v)

def rootKEx(do, n):
	print('calculating rootKEx for n=' + str(n))
	sys.stdout.flush()
	return np.sqrt(do.get3d('KEx', n))

def rootKEy(do, n):
	print('calculating rootKEy for n=' + str(n))
	sys.stdout.flush()
	return np.sqrt(do.get3d('KEy', n))

def rootKEz(do, n):
	print('calculating rootKEz for n=' + str(n))
	sys.stdout.flush()
	return np.sqrt(do.get3d('KEz', n))

def dKE(do, n):
	print('calculating dKE for n=' + str(n))
	sys.stdout.flush()
	dvx = do.get3d('dvx', n)
	dvy = do.get3d('dvy', n)
	dvz = do.get3d('dvz', n)
	rho = do.get3d('rho', n)
	dKE = 0.5*rho*(np.square(dvx) + np.square(dvy) + np.square(dvz))
	return dKE

def rootRhoVx(do, n):
	print('calculating rootRhoVx for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('vx', n)
	rho = do.get3d('rho', n)
	return np.sqrt(0.5*rho)*v

def rootRhoVy(do, n):
	print('calculating rootRhoVy for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('vy', n)
	rho = do.get3d('rho', n)
	return np.sqrt(0.5*rho)*v

def rootRhoVz(do, n):
	print('calculating rootRhoVz for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('vz', n)
	rho = do.get3d('rho', n)
	return np.sqrt(0.5*rho)*v

def rootRhoDvx(do, n):
	print('calculating rootRhoDvx for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('dvx', n)
	rho = do.get3d('rho', n)
	return np.sqrt(0.5*rho)*v

def rootRhoDvy(do, n):
	print('calculating rootRhoDvy for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('dvy', n)
	rho = do.get3d('rho', n)
	return np.sqrt(0.5*rho)*v

def rootRhoDvz(do, n):
	print('calculating rootRhoDvz for n=' + str(n))
	sys.stdout.flush()
	v   = do.get3d('dvz', n)
	rho = do.get3d('rho', n)
	return np.sqrt(0.5*rho)*v


#########################################################################
# FFT STUFF
#########################################################################

def calcPs(do, key, n):
	data       = do.get3d(key, n)
	freqs      = np.fft.fftfreq(data.shape[0], d=do.dx)
	nFreqs     = freqs.shape[0]//2
	freqs      = freqs[:nFreqs]
	fft        = np.fft.fftn(data)[:nFreqs,:nFreqs,:nFreqs]
	ps         = np.square(np.absolute(fft))
	return ps, freqs

def psProfile(do, key, n):
	ps, freqs = calcPs(do, key, n)
	pskx = np.mean(ps, axis=(1,2))
	psky = np.mean(ps, axis=(0,2))
	pskz = np.mean(ps, axis=(0,1))
	psk  = pskx+psky+pskz
	return psk, pskx, psky, pskz, freqs

def psProfileMean(do, key, nStart=None, nEnd=None):
	if nStart is None: nStart = do.nt-10
	if nEnd   is None: nEnd   = do.nt
	pskList  = []
	pskxList = []
	pskyList = []
	pskzList = []
	count    = 0
	for n in range(nStart, nEnd):
		psk, pskx, psky, pskz, freqs = psProfile(do, key, n)
		pskList.append(psk)
		pskxList.append(pskx)
		pskyList.append(psky)
		pskzList.append(pskz)
		count+=1
	psk  = np.asarray(np.mean(pskList,  axis=0))/float(count)
	pskx = np.asarray(np.mean(pskxList, axis=0))/float(count)
	psky = np.asarray(np.mean(pskyList, axis=0))/float(count)
	pskz = np.asarray(np.mean(pskzList, axis=0))/float(count)
	return psk, pskx, psky, pskz, freqs


#########################################################################
# ACF STUFF
#########################################################################

def acf3d(do, key, n):
	data       = do.get3d(key, n)
	nFreqs     = do.nx//2
	fft        = np.fft.fftn(data)
	fftFixed   = np.zeros_like(data)
	ps         = np.square(np.absolute(fft))
	corr       = np.absolute(np.fft.fftn(ps))
	corrNorm   = corr/np.amax(corr)
	return corrNorm

def plotAcf(acfData, figNum=0, extent=[-0.1,0.1,-0.1,0.1]):
	plt.figure(figNum)
	plotData = np.transpose(np.fliplr(acfData))
	plt.imshow(plotData, extent=extent, aspect=1.0, cmap=plt.get_cmap('coolwarm'))
	plt.colorbar()
	plt.clim(-1,1)
	plt.tight_layout()

'''
def acfMean(do, key, shiftMax, ziList, nList, resCut):
	sliceCount = 0
	for n in nList:
		data3d = do.get3d(key, n)
		for zi in ziList:
			print('calculating 2d ACF for'+' n=' + str(n)+' zi='+str(zi))
			sys.stdout.flush()
			data2d = data3d[:,:,zi]
			acf = acf2d(data2d, shiftMax, resCut)
			if sliceCount == 0: acfMean  = acf
			else:               acfMean += acf
			sliceCount+=1
	acfMean /= sliceCount
	return acfMean
'''
'''
def acf4d(do, key, shiftMax, sCut, tCut, resCut):
	nList  = np.arange(0, do.nt, tCut)
	ziList = np.arange(sCut//2, do.nz, sCut)
	acf4d = np.zeros([len(nList), len(ziList), shiftMax*2//resCut+1, shiftMax*2//resCut+1])
	for n in nList:
		data3d = do.get3d(key, n)
		for zi in ziList:
			print('calculating 2d ACF for'+' n=' + str(n)+' zi='+str(zi))
			sys.stdout.flush()
			data2d = data3d[:,:,zi]
			acf    = acf2d(data2d, shiftMax, resCut)
			acf4d[n//tCut, zi//sCut] = acf
	return acf4d
'''
'''
'''
'''
def acf2d(data, shiftMax, resCut):
	corrArray = np.zeros([(2*shiftMax//resCut)+1, (2*shiftMax//resCut)+1])
	array1    = data
	# loop over all desired x and y shifts
	for i in range(-shiftMax//resCut, shiftMax//resCut+1):
		for j in range(-shiftMax//resCut, shiftMax//resCut+1):
			xShift = i * resCut
			yShift = j * resCut
			# the 2d slice of data shifted by xShift and yShift
			array2 = np.roll(np.roll(array1, xShift, axis=0), yShift, axis=1)
			# loop over all rows in the tables and add the correlations together
			tot = 0
			for k in range(0, data.shape[1]):
				# correlate one row of the original and shifted arrays
				corr = np.correlate(array1[:,k], array2[:,k])
				tot = tot + corr
			# assign the final sum to it's place in the corrArray based on the x and y shifts
			corrArray[i+shiftMax//resCut,j+shiftMax//resCut] = tot
	# normalize
	corrArray = corrArray/np.amax(corrArray)
	return corrArray
'''
'''
def acfMean(do, key, shiftMax, ziList, nList, resCut):
	sliceCount = 0
	for n in nList:
		data3d = do.get3d(key, n)
		for zi in ziList:
			print('calculating 2d ACF for'+' n=' + str(n)+' zi='+str(zi))
			sys.stdout.flush()
			data2d = data3d[:,:,zi]
			acf = acf2d(data2d, shiftMax, resCut)
			if sliceCount == 0: acfMean  = acf
			else:               acfMean += acf
			sliceCount+=1
	acfMean /= sliceCount
	return acfMean
'''
'''
def acf4d(do, key, shiftMax, sCut, tCut, resCut):
	nList  = np.arange(0, do.nt, tCut)
	ziList = np.arange(sCut//2, do.nz, sCut)
	acf4d = np.zeros([len(nList), len(ziList), shiftMax*2//resCut+1, shiftMax*2//resCut+1])
	for n in nList:
		data3d = do.get3d(key, n)
		for zi in ziList:
			print('calculating 2d ACF for'+' n=' + str(n)+' zi='+str(zi))
			sys.stdout.flush()
			data2d = data3d[:,:,zi]
			acf    = acf2d(data2d, shiftMax, resCut)
			acf4d[n//tCut, zi//sCut] = acf
	return acf4d
'''
'''
def plotAcf(acfData, figNum=0, extent=[-0.1,0.1,-0.1,0.1]):
	plt.figure(figNum)
	plotData = np.transpose(np.fliplr(acfData))
	plt.imshow(plotData, extent=extent, aspect=1.0, cmap=plt.get_cmap('coolwarm'))
	plt.colorbar()
	plt.clim(-1,1)
	plt.tight_layout()
'''

































#
