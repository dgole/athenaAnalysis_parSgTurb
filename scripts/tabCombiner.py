#!/usr/bin/python
###############################################################################
import numpy as np
import os
import sys
sys.path.append('../python')
import athenaTools as tools
import resource
import time
import multiprocessing
###############################################################################
def getFileNames(basename, timeStep, npc):
	names = [("id"+str(i)+"/"+basename+"-id" + str(i) + "." + tools.getTimeStepString(timeStep) + ".tab") for i in range(npc)]
	names[0] = "id0/"+basename+"."+tools.getTimeStepString(timeStep)+".tab"
	return names
def getFiles(basename, path, names):
	files = [np.loadtxt(path+name) for name in names]
	return files
def combTabs(basename, path, npc, timeStepStart, timeStepEnd):
	for timeStep in range(timeStepStart, timeStepEnd):
		jobName = multiprocessing.current_process().name
		print('process number is ' + jobName)
		print(path + " time step " + str(timeStep))
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
		np.save(outDir+basename+"."+tools.getTimeStepString(timeStep)+".npy", masterArray)
		print("total MB of memory used: " + str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000.0))
		del files
		del masterArray
		del coordsList
		del coordsArray
		del coordsListArray
def checkAndDoLoop(path):
	while True:
		for path in pathList:
			basename="Par_Strat3d"
			outDir=path+'3d/'
			if not os.path.exists(outDir): os.makedirs(outDir)
			contentsList = os.listdir(path)
			npc=0
			for item in contentsList:
				if item[:2]=='id': npc+=1
			contentsList = os.listdir(path+'id0/')
			ntAvail = 0
			for item in contentsList:
				if item[-4:]=='.tab': ntAvail+=1
			contentsList = os.listdir(outDir)
			ntDone = 0
			for item in contentsList:
				if item[-4:]=='.npy': ntDone+=1
			timeStepStart = ntDone
			timeStepEnd   = ntAvail
			print('process ' + multiprocessing.current_process().name + ' checking ' + path)
			print(str(ntAvail) + ' time steps available')
			print(str(ntDone) + ' time steps already done')
			print('combining tabs for indicies ' + str(timeStepStart) + ' through ' + str(timeStepEnd-1))
			combTabs(basename, path, npc, timeStepStart, timeStepEnd)
		time.sleep(10)
###############################################################################
pathList = []
for path in sys.argv[1:]: pathList.append(path)
###############################################################################
print(pathList)
for jobNum in range(len(pathList)):
	print('starting process ' + str(jobNum) + ' for path ' + pathList[jobNum])
	p = multiprocessing.Process(name=str(jobNum), target=checkAndDoLoop, args=(pathList[jobNum],))
	p.start()
























#
