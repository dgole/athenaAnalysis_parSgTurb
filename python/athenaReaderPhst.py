#!/usr/bin/python
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors

####################################################################
# Data class ###############################################
####################################################################
class DataPhst:
	def __init__(self, path, baseName="Par_Strat3d", dt=0.1):
		print("initializing 1d data structure from " + path)
		self.path = path
		fileName  = baseName + '.phst'
		inFile    = open(path+fileName, 'r')
		list1 = []; list2 = []
		for line in inFile:
			if line[0] != '#':
				split = line.split()
				arr   = np.zeros(len(split))
				n = 0
				for s in split:	arr[n] = float(s); n+=1;
				if len(arr) == 11: list1.append(arr)
				if len(arr) == 12: list2.append(arr)
		data1 = np.asarray(list1)
		data2 = np.asarray(list2)
		self.data    = {'t'      : data1[:,0],
						'xavg'   : data2[:,0],
						'yavg'   : data2[:,1],
						'zavg'   : data2[:,2],
						'vxavg'  : data2[:,3],
						'vyavg'  : data2[:,4],
						'vzavg'  : data2[:,5],
						'xvar'   : data2[:,6],
						'yvar'   : data2[:,7],
						'zvar'   : data2[:,8],
						'vxvar'  : data2[:,9],
						'vyvar'  : data2[:,10],
						'vzvar'  : data2[:,11]}
		del list1, list2, arr, split, data1, data2
		#self.header  = {'rho'     : r"$\rho$",
		#				'reynolds': r"$\rho v_x \delta v_y$"}
		self.t      = self.data['t']
		print("phst data of length " + str(self.data['xavg'].shape) + " imported")
	def gettindex(self, t):
		return (np.abs(self.t-t)).argmin()
	#def addCol(self, funcName, key, headerLabel, *args, **kwargs):
	#	print(self.path + ": adding data with key " + key)
	#	self.data[key]   = funcName(self)
	#	self.header[key] = headerLabel

####################################################################
# plotting functions ###############################################
####################################################################
'''
def stPlot(do, key, figNum=0):
	print (do.path + ": making ST plot for key " + key)
	plt.figure(figNum)
	title = do.header[key]
	extent = [0,do.tmax,-do.zmax,do.zmax]
	aspect  = 0.2*do.tmax/do.zmax
	plotData = np.transpose(np.fliplr(do.data[key]))
	if np.amin(plotData) < 0.0:
		cmapType = 'coolwarm'
		maxVal   = np.amax(np.absolute(plotData))
		norm     = colors.SymLogNorm(maxVal/100.0, linscale=2.0)
	else:
		cmapType = 'viridis'
		norm     = colors.LogNorm()
	plt.imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType), norm=norm)
	plt.title(title);
	plt.xlabel(r"$t \Omega$");
	plt.ylabel(r"$z/H$");
	plt.colorbar(shrink=0.5)
	plt.tight_layout()
'''
'''
def profile(do, key, figNum=0, tStart=None, tEnd=None, legendLabel=None):
	print(do.path + ": making profile plot for key " + key)
	plt.figure(figNum)
	if tStart == None: tStart = do.t[-1]/2.0
	if tEnd   == None: tEnd   = do.t[-1]
	nStart   = do.gettindex(tStart)
	nEnd     = do.gettindex(tEnd)
	plotData = np.mean(do.data[key][nStart:nEnd,:], axis=0)
	title    = do.header[key]
	plt.semilogy(do.z, np.absolute(plotData), label=legendLabel)
	plt.ylabel(do.header[key]);
	plt.xlabel(r"$z/H$");
	plt.tight_layout()
'''

def timeEvo(do, key, figNum=0, legendLabel=None):
	print(do.path + ": making timeEvo plot for key " + key)
	plt.figure(figNum)
	plotData = do.data[key]
	if 10.0*np.amin(np.absolute(plotData[2:])) < np.amax(np.absolute(plotData[2:])):
		plt.semilogy(do.t, np.absolute(plotData), label=legendLabel)
	else:
		plt.plot(do.t, np.absolute(plotData), label=legendLabel)
	plt.ylabel(key);
	plt.xlabel(r"$t \Omega$");
	plt.tight_layout()
'''
def dv(do):
	return np.sqrt(2.0*do.data['KE']/do.data['rho'])
'''
