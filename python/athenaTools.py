#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors
import resource

#m.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def getTimeStepString(i):
	if i > 999:	 zstring = ""
	elif i > 99: zstring = "0"
	elif i > 9:  zstring = "00"
	elif i > -1: zstring = "000"
	return zstring+str(i)

def saveAndClear(plotName, figNum=0, dpi=600, bboxOption=1):
	plt.figure(figNum)
	if bboxOption==1: plt.savefig(plotName, bbox_inches='tight', dpi=dpi)
	else:             plt.savefig(plotName, dpi=dpi)
	plt.clf()

def getColor(n, nStart, nEnd):
	span = nEnd-nStart
	nn   = n-nStart
	r = 1.0 - (nn/span)
	g = 0.0
	b = 0.0 + (nn/span)
	return (r,g,b)

def printMem(marker=None):
	if marker is not None: print(marker)
	print("total MB of memory used: " +
		  str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000.0))




def getFileNames(baseName, n, npc):
	names = [("id"+str(i)+"/"+baseName+"-id" + str(i) + "."
			+ getTimeStepString(n) + ".tab") for i in range(npc)]
	names[0] = "id0/"+baseName+"."+getTimeStepString(n)+".tab"
	return names
def getFile(path, name, resultList):
	print('reading in ' + path + name)
	sys.stdout.flush()
	file = np.loadtxt(path+name)
	resultList.append(file)
def assignToArray(file, fileNum, resultList):
	print('assigning file ' + str(fileNum) + ' of ' + str(npc))
	sys.stdout.flush()
	masterArray = np.zeros([coordsListArray[0].shape[0],
							coordsListArray[1].shape[0],
							coordsListArray[2].shape[0],
							cols-3])
	for i in range(file.shape[0]):
		indicies=[(np.abs(coordsListArray[j]-file[i,3+j])).argmin() for j in range(3)]
		masterArray[indicies[0],indicies[1],indicies[2]]=file[i,3:cols]
	resultList.append(masterArray)
def getNpc(path):
	contentsList = os.listdir(path)
	npc=0
	for item in contentsList:
		if item[:2]=='id': npc+=1
	return npc
def getnAvail(path):
	contentsList = os.listdir(path+'id0/')
	ntAvail = 0
	for item in contentsList:
		if item[-4:]=='.tab': ntAvail+=1
	return ntAvail
def getnDone(path):
	contentsList = os.listdir(path+'3d/')
	ntDone = 0
	for item in contentsList:
		if item[-4:]=='.npy': ntDone+=1
	return ntDone
