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
################################################################################
nProcNow=0
nProcTot=0
################################################################################
def getFileNames(baseName, n, npc):
	names = [("id"+str(i)+"/"+baseName+"-id" + str(i) + "." + tools.getTimeStepString(n) + ".tab") for i in range(npc)]
	names[0] = "id0/"+baseName+"."+tools.getTimeStepString(n)+".tab"
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
	contentsList = os.listdir(outDir)
	ntDone = 0
	for item in contentsList:
		if item[-4:]=='.npy': ntDone+=1
	return ntDone
################################################################################

################################################################################
myNPC = int(sys.argv[1])
path  = str(sys.argv[2])
baseName="Par_Strat3d"
outDir = path+'3d/'
if not os.path.exists(outDir): os.makedirs(outDir)
npc    = getNpc(path)
################################################################################
while True:
	npc    = getNpc(path)
	nAvail = getnAvail(path)
	nDone  = getnDone(path)
	print('checking ' + path)
	print(str(nAvail) + ' time steps available')
	print(str(nDone) + ' time steps already done')
	sys.stdout.flush()
	if nAvail > nDone:
		n = nDone
		########################################################################
		# read in files (parallel)
		fileNames = getFileNames(baseName, n, npc)
		manager   = mp.Manager()
		files     = manager.list()
		jobs=[]
		for name in fileNames:
			proc=mp.Process(target=getFile, args=(path,name,files))
			jobs.append(proc)
		nRead=0
		while nRead < npc:
			for n1 in range(myNPC):
				jobs[nRead+n1].start()
			for n1 in range(myNPC):
				jobs[nRead+n1].join()
			nRead+=myNPC
		########################################################################
		# set-up (serial)
		res             = abs(1.0/(files[0][0,3]-files[0][1,3]))
		cols            = files[0].shape[1]
		coordsList      = [file[:,3:6] for file in files]
		coordsArray     = np.asarray(coordsList)
		coordsListArray = [np.unique(coordsArray[:,:,i]) for i in range(3)]
		masterArray     = np.zeros([coordsListArray[0].shape[0],
									coordsListArray[1].shape[0],
									coordsListArray[2].shape[0],
									cols-3])
		masterLength    = len(files)*files[0].shape[0];
		########################################################################
		# assign data to 3d array
		manager = mp.Manager()
		resultList = manager.list()
		jobs=[]; fileNum=0
		for file in files:
			proc=mp.Process(target=assignToArray, args=(file,fileNum,resultList))
			jobs.append(proc)
			fileNum+=1
		nRead=0
		while nRead < npc:
			resultList = manager.list()
			for n1 in range(myNPC):
				jobs[nRead+n1].start()
			for n1 in range(myNPC):
				jobs[nRead+n1].join()
			for result in resultList:
				masterArray+=result
			nRead+=myNPC
		########################################################################
		print('writing master file')
		sys.stdout.flush()
		np.save(outDir+baseName+"."+tools.getTimeStepString(n)+".npy", masterArray)
		print('writing is done')
		sys.stdout.flush()
		del files, masterArray, coordsList, coordsArray, coordsListArray
		time.sleep(5)
	else:
		print('all avaliable output is done')
		sys.stdout.flush()
		break
























#
