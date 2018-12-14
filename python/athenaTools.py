#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
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

def saveAndClear(plotName, figNum=0):
	plt.figure(figNum)
	plt.savefig(plotName, bbox_inches='tight', dpi=600)
	plt.clf()

def getColor(n, nStart, nEnd):
	span = nEnd-nStart
	nn   = n-nStart
	r = 1.0 - (nn/span)
	g = 0.0
	b = 0.0 + (nn/span)
	return (r,g,b)
