#!/usr/bin/python
from __future__ import unicode_literals
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import os
import math
from scipy import fftpack as fft
import sys
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


headerTex1d = [r"$z$", r"$\rho$", r"$P$", r"$KE_x$", r"$KE_y$", r"$KE_z$", r"$KE$", r"$\mathbf{\rho v_x \delta v_y}$", r"$ME_x$", r"$ME_y$", r"$ME_z$", r"$ME$", r"$B_x$", r"$B_y$", r"$B_z$", r"$B_x B_y$"]

class Data1d:
	def __init__(self, path, iStart, iEnd, baseName="SB", dt=0.1):
		print "initializing 1d data structure from " + path
		self.path=path
		self.names = [baseName + "." + getTimeStepString(i) + ".1d" for i in range(iStart, iEnd)]
		self.files = [np.loadtxt(path+self.names[i]) for i in range(iStart, iEnd)]
		self.dataArray = np.asarray(self.files)
		self.data = [self.dataArray[...,...,i] for i in range(self.dataArray.shape[2])]
		del self.files
		del self.dataArray	
		self.nt = self.data[0].shape[0]
		self.nz = self.data[0].shape[1]  
		self.zArray = self.data[0][0,...]
		self.zmax = np.round(-self.data[0][0,0])
		self.tmax = dt*self.nt
		self.tArray = np.arange(0, self.tmax, dt)
		self.pdfName = PdfPages(path + "/analysis1d.pdf")
		self.header = headerTex1d
		print "data of shape " + str(self.data[0].shape) + " imported for " + str(self.nt) + " time steps"
		print "total MB of memory used: " + str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000.0)

	def getzindex(self, z):
		return (np.abs(self.zArray-z)).argmin()

	def gettindex(self, t):
		return (np.abs(self.tArray-t)).argmin()

	def addCol(self, funcName, label, *args, **kwargs):
		print self.path + ": adding column named " + label
		self.data.append(funcName(self.data))
		self.header.append(label)		

	def stPlot(self, col, cmapType="viridis", logOption=0, save=None, savePath=None, clim1=None, clim2=None):
		print self.path + ": making ST plot for column " + str(col)
		if savePath is None:
			savePath=self.path
		if logOption==0:
			plotData = self.data[col]; title = self.header[col];
		if logOption==1:
			plotData = np.log10(np.absolute(self.data[col])); title = "log " + self.header[col];
		plt.imshow(np.transpose(np.fliplr(plotData)), extent=[0,self.tmax,-self.zmax,self.zmax], aspect=(0.2*self.tmax/self.zmax), cmap=plt.get_cmap(cmapType))
		plt.title(title); plt.xlabel("Time (orbits)"); plt.ylabel(r"$z/H$");
		plt.colorbar(shrink=0.5)
		if (clim1 is not None and clim2 is not None):
			plt.clim(clim1,clim2)
		if (clim1 is not None and clim2 is None):
			plt.clim(clim1,np.amax(plotData))
		if (clim1 is None and clim2 is not None):
			plt.clim(np.aminx(plotData),clim2)
		plt.tight_layout()
		if save=="png":
			plt.savefig(savePath + "ST_" + str(col) + ".png", bbox_inches='tight')
			print "saved ST plot for column " + str(col) + " to png"
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print "saved ST plot for column " + str(col) + " to pdf"
		plt.clf()
		
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


def vx(data):
	vx = np.power(2.0*data[3]/data[1],0.5)
	return vx

def vy(data):
	vy = np.power(2.0*data[4]/data[1],0.5)
	return vy

def vz(data):
	vz = np.power(2.0*data[5]/data[1],0.5)
	return vz






