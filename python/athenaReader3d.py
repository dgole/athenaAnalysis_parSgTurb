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
		self.xmax   = np.round(-self.x[0],1)
		self.ymax   = np.round(-self.y[0],1)
		self.zmax   = np.round(-self.z[0],1)
		self.dt     = dt
		self.tmax   = dt*self.nt
		self.t      = np.arange(0, self.tmax-dt/2.0, dt)
		print(self.t)
		self.q      = 1.5
		self.omega  = 1.0
		self.vk     = -self.q*self.omega*self.x
		#self.shearTime = (self.y[0]-self.y[-1])/self.vk[-1]
		#print(self.shearTime)
		self.shearTime = 1/(self.q*self.omega)
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
	def get3d_unsheared(self, key, n, *args, **kwargs):
		data       = self.get3d(key, n)
		newData    = np.zeros_like(data)
		tSinceEven = np.mod(self.t[n], self.shearTime)
		shearDisp  = (tSinceEven/self.shearTime) * self.x * 2.0
		#print(shearDisp)
		shearNumFloat = np.floor((shearDisp/self.dx)+0.5)
		#print(shearNumFloat)
		for i in range(data.shape[0]):
			newData[i,:,:] = np.roll(data[i,:,:], -int(shearNumFloat[i]), axis=0)
			#print(i, n, -shearI)
		return newData
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

def profile(do, key, figNum=0, tStart=None, tEnd=None, legendLabel=None, absAvg=1, absPlot=1, color='b'):
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
		else:         data3d =             do.get3d(key, n)
		data1d    = np.mean(data3d, axis=(0,1))
		plotData += data1d
		nCount   += 1
	plotData /= nCount
	title = do.header[key]
	if absPlot==1:
		plt.plot(do.z, np.absolute(plotData), label=legendLabel, color=color)
	else:
		plt.plot(do.z, plotData , label=legendLabel, color=color)
	plt.ylabel(do.header[key]);
	plt.xlabel(r"$z/H$");
	plt.tight_layout()

def timeEvo(do, key, figNum=0, legendLabel=None, absAvg=1, absPlot=1, color='b'):
	print(do.path + ": making timeEvo plot for key " + key)
	sys.stdout.flush()
	plt.figure(figNum)
	plotData = np.zeros(do.nt)
	for n in range(0, do.nt):
		if absAvg==1: plotData[n] = np.mean(np.absolute(do.get3d(key, n)))
		else:         plotData[n] = np.mean(            do.get3d(key, n))
	plotData /= do.nt
	title = do.header[key]
	if absPlot==1:
		plt.semilogy(do.t, np.absolute(plotData), label=legendLabel, color=color)
	else:
		plt.plot(do.t, plotData , label=legendLabel, color=color)
	plt.ylabel(do.header[key]);
	plt.xlabel(r"$t \Omega$");
	plt.tight_layout()

def slicePlot(do, key, n=None, figNum=0, axis='z', coord=0.0):
	print(do.path + ": making slice plot for key " + key)
	if n==None: n = do.nt - 10
	data3d = do.get3d(key, n)
	if axis=='x':
		index  = do.getxindex(coord)
		data2d = data3d[index, :, :]
		extent = [-do.ymax, do.ymax, -do.zmax, do.zmax]
		xlabel = r"$y/H$"
		ylabel = r"$z/H$"
	if axis=='y':
		index  = do.getyindex(coord)
		data2d = data3d[:, index, :]
		extent = [-do.xmax, do.xmax, -do.zmax, do.zmax]
		xlabel = r"$x/H$"
		ylabel = r"$z/H$"
	if axis=='z':
		index  = do.getzindex(coord)
		data2d = data3d[:, :, index]
		extent = [-do.xmax, do.xmax, -do.ymax, do.ymax]
		xlabel = r"$x/H$"
		ylabel = r"$y/H$"
	title = do.header[key]
	aspect  = 1.0
	plotData = np.transpose(np.fliplr(data2d))
	if np.amin(plotData) < 0.0:
		cmapType = 'coolwarm'
		maxVal   = np.amax(np.absolute(plotData))
		#norm     = colors.SymLogNorm(maxVal/100.0, linscale=2.0)
	else:
		cmapType = 'viridis'
		#norm     = colors.LogNorm()
	#plt.imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType), norm=norm)
	plt.imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.colorbar()
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
	aAvg1d = np.mean(np.absolute(a), axis=(0,1))
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
# PSPEC STUFF
#########################################################################

def calcPs(do, key, n, unsheared=True):
	if unsheared: data = do.get3d_unsheared(key, n)
	else:         data = do.get3d(key, n)
	freqs      = np.fft.fftfreq(data.shape[0], d=do.dx)
	nFreqs     = freqs.shape[0]//2
	freqs      = freqs[:nFreqs]
	fft        = np.fft.fftn(data)[:nFreqs,:nFreqs,:nFreqs]
	ps         = np.square(np.absolute(fft))
	return ps, freqs

def psProfiles(do, key, n, unsheared=True):
	ps, freqs = calcPs(do, key, n, unsheared=unsheared)
	psk       = np.zeros(ps.shape[0])
	count     = np.zeros_like(psk)
	for i in range(ps.shape[0]):
		for j in range(ps.shape[1]):
			for k in range(ps.shape[2]):
				dist          = np.sqrt(i*i+j*j+k*k)
				index         = int(np.floor(dist))
				if index < freqs.shape[0]:
					psk[index]   += ps[i,j,k]
					count[index] += 1
	psk /= count
	return psk, freqs

def psProfileMean(do, key, nStart=None, nEnd=None, unsheared=True):
	#if nStart is None: nStart = do.nt // 2
	if nStart is None: nStart = do.nt // 2
	if nEnd   is None: nEnd   = do.nt
	pskList = []
	count   = 0
	for n in range(nStart, nEnd):
		if np.mod(do.t[n], do.shearTime)<1.e-8:
			psk, freqs = psProfiles(do, key, n, unsheared=unsheared)
			pskList.append(psk)
			count += 1
	psk = np.asarray(np.mean(pskList,  axis=0))/float(count)
	return psk, freqs


#########################################################################
# ACF STUFF
#########################################################################

def acf3d(do, key, n, unsheared=True):
	print('calculating ACF for'+' n=' + str(n))
	sys.stdout.flush()
	if unsheared: data = do.get3d_unsheared(key, n)
	else:         data = do.get3d(key, n)
	nFreqs     = do.nx//2
	fft        = np.fft.fftn(data)
	fftFixed   = np.zeros_like(data)
	ps         = np.square(np.absolute(fft))
	corr       = np.absolute(np.fft.fftn(ps))
	corrNorm   = corr/np.amax(corr)
	return corrNorm

def acf3dto2d_xy(do, key, n, unsheared=True):
	nn        = do.nx//2
	acf       = acf3d(do, key, n, unsheared=unsheared)
	acf2d     = np.mean(acf, axis=2)
	acf2d2   = np.zeros_like(acf2d)
	acf2d2[:nn,:nn] = acf2d[nn:,nn:] # 3 into 1
	acf2d2[:nn,nn:] = acf2d[nn:,:nn] # 4 into 2
	acf2d2[nn:,nn:] = acf2d[:nn,:nn] # 1 into 3
	acf2d2[nn:,:nn] = acf2d[:nn,nn:] # 2 into 4
	return acf2d2

def acfMean(do, key, nStart=None, nEnd=None, unsheared=True):
	if nStart is None: nStart = do.nt // 2
	if nEnd   is None: nEnd   = do.nt
	for n in range(nStart,nEnd):
		acf2d = acf3dto2d_xy(do, key, n, unsheared=unsheared)
		if n==nStart: acfMean  = acf2d
		else:         acfMean += acf2d
	acfMean /= np.amax(acfMean)
	return acfMean

def plotAcf(acf, figNum=0, extent=[-0.1,0.1,-0.1,0.1]):
	plt.figure(figNum)
	plotData = np.transpose(np.fliplr(acf))
	plt.imshow(plotData, extent=extent, aspect=1.0, cmap=plt.get_cmap('viridis'))
	plt.clim(0,1)
	plt.colorbar()
	plt.tight_layout()

def plotAcfDiverging(acf, figNum=0, extent=[-0.1,0.1,-0.1,0.1]):
	plt.figure(figNum)
	plotData = np.transpose(np.fliplr(acf))
	plt.imshow(plotData, extent=extent, aspect=1.0, cmap=plt.get_cmap('coolwarm'))
	plt.clim(-1,1)
	plt.colorbar()
	plt.tight_layout()
































#
