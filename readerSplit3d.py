#!/usr/bin/python
from __future__ import unicode_literals
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
import resource
import os
m.rcParams['text.usetex'] = True
m.rcParams['text.latex.unicode'] = True

def getTimeStepString(i):
	if i > 999:	 zstring = ""
	elif i > 99: zstring = "0"
	elif i > 9:  zstring = "00"
	elif i > -1: zstring = "000"
	return zstring+str(i)

class Data3d:
	def __init__(self, path, runID, checkKey='By', baseName="StratCooling", dt=1.0):
		print("initializing 3d data structure from " + path)
		self.runID         = runID
		self.path          = path
		self.baseName      = baseName
		self.dt            = dt
		self.baseNames     = [self.path + self.baseName + "." + getTimeStepString(n) + '.' for n in range(0, 10000)]
		self.nt            = 0
		for n in range(0,10000):
			if os.path.isfile(self.baseNames[n] + checkKey + '.npy'): self.nt += 1;
			else: break;
		self.t             = np.arange(0, (self.nt+1)*dt, dt)
		self.x             = np.load(self.baseNames[0] + 'x' + '.npy');
		self.y             = np.load(self.baseNames[0] + 'y' + '.npy');
		self.z             = np.load(self.baseNames[0] + 'z' + '.npy'); 
		self.nx            = (self.x.shape)[0]
		self.ny            = (self.y.shape)[0]
		self.nz            = (self.z.shape)[0]
		self.header        = {'rho' : r"$\rho$",
						              'vx'  : r"$v_x$" ,
						              'vy'  : r"$v_y$" ,
						              'vz'  : r"$v_z$" ,
						              'P'   : r"$P$"   ,
						              'Bx'  : r"$B_x$" ,
						              'By'  : r"$B_y$" ,
						              'Bz'  : r"$B_z$"  }
	def getData3d(self, key, n):
		fileName = self.baseNames[n] + key + '.npy'
		print('reading ' + fileName)
		returnArray = np.load(fileName)
		print("total MB of memory used: " + str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000.0))
		return returnArray


















