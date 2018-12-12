#!/usr/bin/python
###############################################################################
import numpy as np
import os
import sys
sys.path.append('../python')
import athenaTools as tools
import resource
import time
import multiprocessing as mp
nProc=0
nDoingDict=dict()
###############################################################################
def getFileNames(basename, timeStep, npc):
	names = [("id"+str(i)+"/"+basename+"-id" + str(i) + "." + tools.getTimeStepString(timeStep) + ".tab") for i in range(npc)]
	names[0] = "id0/"+basename+"."+tools.getTimeStepString(timeStep)+".tab"
	return names
def getFiles(basename, path, names):
	print("reading in files")
	sys.stdout.flush()
	files = []
	fileNum = 0
	for name in names:
		print('reading in file number ' + str(fileNum) + ' of ' + str(len(names)))
		print("total MB of memory used: " + str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000.0))
		files.append(np.loadtxt(path+name))
		fileNum+=1
	return files
def combTabs(basename, path, npc, timeStep):
	global nProc
	global nDoingDict
	outDir = path+'3d/'
	# read in files
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
	fileNum = 0
	for file in files:
		print('shuffling file number ' + str(fileNum) + ' of ' + str(len(files)))
		print("total MB of memory used: " + str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000.0))
		sys.stdout.flush()
		fileNum+=1
		for i in range(file.shape[0]):
			indicies=[(np.abs(coordsListArray[j]-file[i,3+j])).argmin() for j in range(3)]
			masterArray[indicies[0],indicies[1],indicies[2]]=file[i,3:cols]
	# export table
	np.save(outDir+basename+"."+tools.getTimeStepString(timeStep)+".npy", masterArray)
	del files
	del masterArray
	del coordsList
	del coordsArray
	del coordsListArray
	nProc-=1
	nDoingDict[path]-=1
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
	contentsList = os.listdir(outDir)
	ntDone = 0
	for item in contentsList:
		if item[-4:]=='.npy': ntDone+=1
	return ntDone

###############################################################################
myNPC = int(sys.argv[1])
basename="Par_Strat3d"
pathList = []
for path in sys.argv[2:]:
	pathList.append(path)
	nDoingDict[path]=0
###############################################################################
while True:
	for path in pathList:
		nDoing = nDoingDict[path]
		outDir = path+'3d/'
		if not os.path.exists(outDir): os.makedirs(outDir)
		npc     = getNpc(path)
		nAvail = getnAvail(path)
		nDone  = getnDone(path)
		print(' ')
		print('######################################################################')
		print(str(nProc) + ' processes currently running')
		print('checking ' + path)
		print(str(nAvail) + ' time steps available')
		print(str(nDone) + ' time steps already done')
		print(str(nDoing) + ' time steps currently being done')
		sys.stdout.flush()
		if nAvail > nDone + nDoing:
			timeStep = nDone + nDoing
			if myNPC>1:
				while True:
					if nProc < myNPC:
						print('starting process to do ' + path + ', n=' + str(timeStep))
						sys.stdout.flush()
						p = mp.Process(target=combTabs, args=(basename, path, npc, timeStep))
						nProc +=1
						nDoingDict[path] +=1
						break;
					else:
						print('at max number of processes, waiting...')
						sys.stdout.flush()
						time.sleep(30)
			else:
				print('starting ' + path + ', n=' + str(timeStep))
				sys.stdout.flush()
				combTabs(basename, path, npc, timeStep)

	time.sleep(2)
























#
