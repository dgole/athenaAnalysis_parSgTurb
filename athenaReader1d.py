#!/usr/bin/python
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors

#m.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def getTimeStepString(i):
	if i > 999:	 zstring = ""
	elif i > 99: zstring = "0"
	elif i > 9:  zstring = "00"
	elif i > -1: zstring = "000"
	return zstring+str(i)

#headerTex1d = [r"$z$", r"$\rho$", r"$P$", r"$T$", r"$E_x$", r"$E_tot$", r"$KE_x$", r"$KE_y$", r"$KE_z$", r"$KE$", r"$\mathbf{\rho v_x \delta v_y}$", r"$ME_x$", r"$ME_y$", r"$ME_z$", r"$ME$", r"$B_x$", r"$B_y$", r"$B_z$", r"$B_x B_y$"]
headerTex1d = [r"$z$", r"$\rho$", r"$P$", r"$KE_x$", r"$KE_y$", r"$KE_z$", r"$KE$", r"$\rho v_x \delta v_y$"]

class Data:
	def __init__(self, path, baseName="Par_Strat3d", dt=0.1):
		print("initializing 1d data structure from " + path)
		self.path      = path
		fileList       = os.listdir(self.path)
		self.nFiles    = len(fileList)
		self.names     = [baseName + "." + getTimeStepString(n) + ".1d" for n in range(self.nFiles)]
		self.files     = [np.loadtxt(path+self.names[n]) for n in range(self.nFiles)]
		self.dataArray = np.asarray(self.files)
		self.data      = [self.dataArray[:,:,i] for i in range(self.dataArray.shape[2])]
		del self.files
		del self.dataArray
		self.nt     = self.data[0].shape[0]
		self.nz     = self.data[0].shape[1]
		self.z      = self.data[0][0,:]
		self.zmax   = np.round(-self.data[0][0,0],1)
		self.tmax   = dt*self.nt
		self.t      = np.arange(0, self.tmax, dt)
		self.header = headerTex1d
		print("data of shape " + str(self.data[0].shape) + " imported")

	def getzindex(self, z):
		return (np.abs(self.z-z)).argmin()

	def gettindex(self, t):
		return (np.abs(self.t-t)).argmin()

	#def addCol(self, funcName, label, *args, **kwargs):
		#print self.path + ": adding column named " + label
		#self.data.append(funcName(self.data))
		#self.header.append(label)


def stPlot(do, col, figNum=0, cmapType="viridis", logOption=0):
	print (do.path + ": making ST plot for column " + str(col))
	plt.figure(figNum)
	title = do.header[col];
	plotData = do.data[col]
	if np.amin(plotData) < 0.0 and cmapType == 'viridis':
		plotData = np.absolute(plotData)
	if cmapType == 'viridis':
		plt.imshow(
				   np.transpose(np.fliplr(plotData)),
				   extent=[0,do.tmax,-do.zmax,do.zmax],
				   aspect=(0.2*do.tmax/do.zmax),
				   cmap=plt.get_cmap(cmapType),
				   norm=colors.LogNorm()
				   )
	elif cmapType == 'coolwarm':
		maxVal = np.amax(np.absolute(plotData))
		plt.imshow(
				   np.transpose(np.fliplr(plotData)),
				   extent=[0,do.tmax,-do.zmax,do.zmax],
				   aspect=(0.2*do.tmax/do.zmax),
				   cmap=plt.get_cmap(cmapType),
				   norm=colors.SymLogNorm(maxVal/100.0, linscale=1.0)
				   )
	plt.title(title);
	plt.xlabel(r"$t \Omega$");
	plt.ylabel(r"$z/H$");
	plt.colorbar(shrink=0.5)
	plt.tight_layout()


'''
def profile(self, col, tStart, tEnd, logOption=0, save=None, savePath=None, legendLabel=None):
	print self.path + ": making profile plot for column " + str(col)
	nStart=self.gettindex(tStart); nEnd=self.gettindex(tEnd);
	if savePath is None:
		savePath=self.path
	plotData = np.mean(self.data[col][nStart:nEnd,...], axis=0); title = self.header[col];
	if logOption==0 and save is not None:
		plt.plot(self.zArray, plotData);	plt.ylabel(self.header[col]);
	if logOption==1 and save is not None:
		plt.semilogy(self.zArray, np.absolute(plotData));	plt.ylabel(self.header[col]);
	if logOption==0 and save is None:
		plt.plot(self.zArray, plotData, label=legendLabel);	plt.ylabel(self.header[col]);
	if logOption==1 and save is None:
		plt.semilogy(self.zArray, np.absolute(plotData), label=legendLabel);	plt.ylabel(self.header[col]);
	plt.xlabel(r"$z/H$");
	plt.tight_layout()
	if save=="png":
		plt.savefig(savePath + "profile_" + str(col) + ".png", bbox_inches='tight')
		print "saved profile plot for column " + str(col) + " to png"
	if save=="pdf":
		plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
		print "saved profile plot for column " + str(col) + " to pdf"
	if save is None:
		#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.legend()
	if save is not None:
		plt.clf()

def timeEvo(self, col, zStart, zEnd, logOption=0, save=None, savePath=None, legendLabel=None):
	print self.path + ": making timeEvo plot for column " + str(col)
	ziStart=self.getzindex(zStart)
	ziEnd=self.getzindex(zEnd)
	if savePath is None:
		savePath=self.path
	plotData = np.mean(self.data[col][...,ziStart:ziEnd], axis=1); title = self.header[col];
	if logOption==0 and save is not None:
		plt.plot(self.tArray, plotData); plt.ylabel(self.header[col]);
	if logOption==1 and save is not None:
		plt.semilogy(self.tArray, np.absolute(plotData)); plt.ylabel(self.header[col]);
	if logOption==2 and save is not None:
		plt.loglog(self.tArray, np.absolute(plotData));	plt.ylabel("log " + self.header[col]);
	if logOption==0 and save is None:
		plt.plot(self.tArray, plotData, label=legendLabel); plt.ylabel(self.header[col]);
	if logOption==1 and save is None:
		plt.semilogy(self.tArray, np.absolute(plotData), label=legendLabel); plt.ylabel(self.header[col]);
	if logOption==2 and save is None:
		plt.loglog(self.tArray, np.absolute(plotData), label=legendLabel);	plt.ylabel("log " + self.header[col]);
	plt.xlabel("Time (orbits)");
	plt.tight_layout()
	if save=="png":
		plt.savefig(savePath + "timeEvo_" + str(col) + ".png", bbox_inches='tight')
		print "saved timeEvo plot for column " + str(col) + " to png"
	if save=="pdf":
		plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
		print "saved timeEvo plot for column " + str(col) + " to pdf"
	if save is None:
		plt.legend()
	if save is not None:
		plt.clf()
'''





























'''
	def makeStandardPlotsPdf(self, t1, t2, t3):
		#density
		col=1
		self.stPlot(col, logOption=1, save="pdf")
		self.profile(col, t2, t3, logOption=1, save="pdf")
		#pressure
		col=2
		self.stPlot(col, logOption=1, save="pdf")
		self.profile(col, t2, t3, logOption=1, save="pdf")
		#temp
		col=3
		self.stPlot(col, clim1=0.4, clim2=1.2, save="pdf")
		self.profile(col, t2, t3, save="pdf")
		#KE
		col=9
		self.stPlot(col, logOption=1, save="pdf", clim1=-7.0)
		self.profile(col, t2, t3, logOption=1, save="pdf")
		self.timeEvo(col, t1, t3, logOption=1, save="pdf")
		#Reynolds
		col=10
		self.stPlot(col, logOption=1, save="pdf", clim1=-7.0)
		self.profile(col, t2, t3, logOption=1, save="pdf")
		self.timeEvo(col, t1, t3, logOption=1, save="pdf")
		#ME
		col=14
		self.stPlot(col, logOption=1, save="pdf", clim1=-7.0)
		self.profile(col, t2, t3, logOption=1, save="pdf")
		self.timeEvo(col, t1, t3, logOption=1, save="pdf")
		#Maxwell
		col=18
		self.stPlot(col, logOption=1, save="pdf", clim1=-7.0)
		self.profile(col, t2, t3, logOption=1, save="pdf")
		self.timeEvo(col, t1, t3, logOption=1, save="pdf")







def dv1d(data):
	return np.sqrt(2.0*data[9]/data[1])/np.sqrt(data[3])

def transEff(data):
	return data[10]/data[9]

def btot3d(data):
	return np.power(np.power(data[8],2)+np.power(data[9],2)+np.power(data[10],2),0.5)

def maxwell3d(data):
	return -data[8]*data[9]

def vx(data):
	return -data[8]*data[9]

def dv3d(data):
	dv=np.zeros_like(data[0]);
	for n in range(0,data[0].shape[0]):
		shiftedArrayListVx = []
		shiftedArrayListVy = []
		shiftedArrayListVz = []
		for iShift in range(-2, 3):
			for jShift in range(-2, 3):
				for kShift in range(-2, 3):
					shiftedArrayListVx.append(np.roll(np.roll(np.roll(data[4][n], iShift, axis=0), jShift, axis=1), kShift, axis=2))
					shiftedArrayListVy.append(np.roll(np.roll(np.roll(data[5][n], iShift, axis=0), jShift, axis=1), kShift, axis=2))
					shiftedArrayListVz.append(np.roll(np.roll(np.roll(data[6][n], iShift, axis=0), jShift, axis=1), kShift, axis=2))
		shiftedArrayArrayVx=np.asarray(shiftedArrayListVx)
		shiftedArrayArrayVy=np.asarray(shiftedArrayListVy)
		shiftedArrayArrayVz=np.asarray(shiftedArrayListVz)
		del shiftedArrayListVx
		del shiftedArrayListVy
		del shiftedArrayListVz
		dv[n] = np.var(shiftedArrayArrayVx, axis=0) + np.var(shiftedArrayArrayVy, axis=0) + np.var(shiftedArrayArrayVz, axis=0)
		del shiftedArrayArrayVx
		del shiftedArrayArrayVy
		del shiftedArrayArrayVz
		print "n=" + str(n) + " total MB of memory used: " + str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000.0)
	return np.sqrt(dv)/np.sqrt(data[7]/data[3])

def drho3d(data):
	drho=np.zeros_like(data[0]);
	for n in range(0,data[0].shape[0]):
		shiftedArrayList = []
		for iShift in range(-2, 3):
			for jShift in range(-2, 3):
				for kShift in range(-2, 3):
					shiftedArrayList.append(np.roll(np.roll(np.roll(data[3][n], iShift, axis=0), jShift, axis=1), kShift, axis=2))
		shiftedArrayArray=np.asarray(shiftedArrayList)
		del shiftedArrayList
		drho[n] = np.var(shiftedArrayArray, axis=0)
		del shiftedArrayArray
		print "n=" + str(n) + " total MB of memory used: " + str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000.0)
	return np.sqrt(drho)/data[3]

def zeros3d(data):
	return np.zeros_like(data[0]);

def keFluxx(data):
	ke = 0.5*data[3]*(np.power(data[4], 2) + np.power(data[5], 2) + np.power(data[6], 2))
	keFlux = data[4]*ke
	return keFlux

def keFluxy(data):
	ke = 0.5*data[3]*(np.power(data[4], 2) + np.power(data[5], 2) + np.power(data[6], 2))
	keFlux = data[5]*ke
	return keFlux

def keFluxz(data):
	ke = 0.5*data[3]*(np.power(data[4], 2) + np.power(data[5], 2) + np.power(data[6], 2))
	keFlux = data[6]*ke
	return keFlux
'''
