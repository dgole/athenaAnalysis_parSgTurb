#!/usr/bin/python
import numpy as np
import matplotlib as m
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
		self.cols    = {'rho'    : 3,
						'vx'     : 4,
						'vy'     : 5,
						'vz'     : 6,
						'dpar'   : 7,
						'vxpar'  : 8,
						'vypar'  : 9,
						'vypar'  : 10}
		self.funcs   = {
						'drho'     : drho,
						'dvx'      : dvx,
						'dvy'      : dvy,
						'dvz'      : dvz,
						'dv'       : dv,
						'drhoNorm' : drhoNorm,
						'dvxNorm'  : dvxNorm,
						'dvyNorm'  : dvyNorm,
						'dvzNorm'  : dvzNorm,
						'dvNorm'   : dvNorm
									         	}
		#self.header  = {'rho'     : r"$\rho$",
		#				'P'       : r"$P$",
		#				'KEx'     : r"$KE_x$",
		#				'KEy'     : r"$KE_y$",
		#				'KEz'     : r"$KE_z$",
		#				'KE'      : r"$KE$",
		#				'reynolds': r"$\rho v_x \delta v_y$"}
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
	daNorm = (a - aAvg3d) / np.absolute(a)
	return daNorm


#########################################################################
# Callable with do.get3d()
#########################################################################

def drho(do, n):
	print('calculating drho for n=' + str(n))
	da = getPert(do.get3d('rho', n))
	return da

def drhoNorm(do, n):
	print('calculating drhoNorm for n=' + str(n))
	daNorm = getPertNorm(do.get3d('rho', n))
	return daNorm

def dvx(do, n):
	print('calculating dvx for n=' + str(n))
	da = getPert(do.get3d('vx', n))
	return da

def dvxNorm(do, n):
	print('calculating dvxNorm for n=' + str(n))
	daNorm = getPertNorm(do.get3d('vx', n))
	return daNorm

def dvy(do, n):
	print('calculating dvy for n=' + str(n))
	da = getPert(do.get3d('vy', n))
	return da

def dvyNorm(do, n):
	print('calculating dvyNorm for n=' + str(n))
	daNorm = getPertNorm(do.get3d('vy', n))
	return daNorm

def dvz(do, n):
	print('calculating dvz for n=' + str(n))
	da = getPert(do.get3d('vz', n))
	return da

def dvzNorm(do, n):
	print('calculating dvzNorm for n=' + str(n))
	daNorm = getPertNorm(do.get3d('vz', n))
	return daNorm

def dv(do, n):
	print('calculating dv for n=' + str(n))
	dvx = do.get3d('dvx', n)
	dvy = do.get3d('dvy', n)
	dvz = do.get3d('dvz', n)
	dv  = np.sqrt(np.square(dvx) + np.square(dvy) + np.square(dvz))
	return dv

def dvNorm(do, n):
	print('calculating dvNorm for n=' + str(n))
	dvx = do.get3d('dvx', n)
	dvy = do.get3d('dvy', n)
	dvz = do.get3d('dvz', n)
	dv  = np.sqrt(np.square(dvx) + np.square(dvy) + np.square(dvz))
	v   = np.sqrt(np.square(do.get3d('vx',n)) + np.square(do.get3d('vy',n)) + np.square(do.get3d('vz',n)))
	dvNorm = dv/v
	return dvNorm


#########################################################################
# ACF STUFF
#########################################################################

def acf2d(data, shiftMax):
	corrArray = np.zeros([2*shiftMax+1, 2*shiftMax+1])
	array1    = data
	# loop over all desired x and y shifts
	for i in range(-shiftMax, shiftMax+1):
		for j in range(-shiftMax, shiftMax+1):
			xShift = i
			yShift = j
			# the 2d slice of data shifted by xShift and yShift
			array2 = np.roll(np.roll(array1, xShift, axis=0), yShift, axis=1)
			# loop over all rows in the tables and add the correlations together
			tot = 0
			for k in range(0, data.shape[1]):
				# correlate one row of the original and shifted arrays
				corr = np.correlate(array1[:,k], array2[:,k])
				tot = tot + corr
			# assign the final sum to it's place in the corrArray based on the x and y shifts
			corrArray[i+shiftMax,j+shiftMax] = tot
	# normalize
	corrArray = corrArray/np.amax(corrArray)
	return corrArray

def acfMean(do, key, shiftMax, ziList, nList):
	sliceCount = 0
	for n in nList:
		data3d = do.get3d(key, n)
		for zi in ziList:
			print('calculating 2d ACF for'+' n=' + str(n)+' zi='+str(zi))
			data2d = data3d[:,:,zi]
			acf = acf2d(data2d, shiftMax)
			if sliceCount == 0: acfMean  = acf
			else:               acfMean += acf
			sliceCount+=1
	acfMean /= sliceCount
	return acfMean

def acf4d(do, key, shiftMax):
	acf4d = np.zeros([do.nt, do.nz, shiftMax*2+1, shiftMax*2+1])
	for n in range(40, do.nt, 1):
		data3d = do.get3d(key, n)
		for zi in range(0, do.nz, 1):
			print('calculating 2d ACF for'+' n=' + str(n)+' zi='+str(zi))
			data2d       = data3d[:,:,zi]
			acf          = acf2d(data2d, shiftMax)
			acf4d[n, zi] = acf
	return acf4d

def plotAcf(acfData, figNum=0):
	plt.figure(figNum)
	plotData = np.transpose(np.fliplr(acfData))
	plt.imshow(plotData, extent=[-1,1,-1,1], aspect=1.0, cmap=plt.get_cmap('coolwarm'))
	plt.colorbar()
	plt.clim(-1,1)
	plt.tight_layout()
