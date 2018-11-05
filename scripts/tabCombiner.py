#!/usr/bin/python

####################################################
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import resource
####################################################

# parameters
npc=int(sys.argv[1])
timeStepStart=int(sys.argv[2])
timeStepEnd=int(sys.argv[3])
path=str(sys.argv[4])
outDir=str(sys.argv[5])
basename="Par_Strat3d"

def getTimeStepString(i):
	if i > 999:	 zstring = ""
	elif i > 99: zstring = "0"
	elif i > 9:  zstring = "00"
	elif i > -1: zstring = "000"
	return zstring+str(i)
def getFileNames(basename, timeStep, npc):
	names = [("id"+str(i)+"/"+basename+"-id" + str(i) + "." + getTimeStepString(timeStep) + ".tab") for i in range(npc)] 
	names[0] = "id0/"+basename+"."+getTimeStepString(timeStep)+".tab"
	return names
def getFiles(basename, path, names):
	files = [np.loadtxt(path+name) for name in names]
	return files

for timeStep in range(timeStepStart, timeStepEnd):
	print("combining tab files for " + path + "  time step = " + str(timeStep))
	# read in files
	print("reading in files")
	sys.stdout.flush()
	files=getFiles(basename, path, getFileNames(basename, timeStep, npc))
	res=abs(1.0/(files[0][0,3]-files[0][1,3]))
	cols=files[0].shape[1]
	coordsList=[file[:,3:6] for file in files]
	coordsArray=np.asarray(coordsList)
	coordsListArray=[np.unique(coordsArray[:,:,i]) for i in range(3)]
	# assign data to 3d array
	masterArray=np.zeros([coordsListArray[0].shape[0],coordsListArray[1].shape[0],coordsListArray[2].shape[0],cols-3])
	masterLength=len(files)*files[0].shape[0];
	print("arranging files into 3d numpy array")
	sys.stdout.flush()
	for file in files:
		for i in range(file.shape[0]): 
			indicies=[(np.abs(coordsListArray[j]-file[i,3+j])).argmin() for j in range(3)]
			masterArray[indicies[0],indicies[1],indicies[2]]=file[i,3:cols]
	# export table
	if not os.path.exists(outDir): os.makedirs(outDir)
	np.save(outDir+basename+"."+getTimeStepString(timeStep)+".npy", masterArray)
	print("total MB of memory used: " + str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000.0))
	del files
	del masterArray
	del coordsList
	del coordsArray
	del coordsListArray

















































