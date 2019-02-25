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
myNPC    = int(sys.argv[1])
path     = str(sys.argv[2])
baseName = "Par_Strat3d"
outDir   = path+'3d/'
if not os.path.exists(outDir): os.makedirs(outDir)

################################################################################

def processFile(path, name, resultList):
	#print('processing ' + path + name)
	sys.stdout.flush()
	file = np.loadtxt(path+name)
	coordsList      = file[:,3:6]
	coordsArray     = np.asarray(coordsList)
	coordsListArray = [np.unique(coordsArray[:,i]) for i in range(3)]
	data1 = file.reshape((len(coordsListArray[0]),
		 				  len(coordsListArray[1]),
						  len(coordsListArray[2]),
						  file.shape[1]),
						  order='F')
	data = data1[:,:,:,3:]
	del file, data1
	resultList.append(data)

################################################################################

while True:
	npc    = tools.getNpc(path)
	nAvail = tools.getnAvail(path)
	nDone  = tools.getnDone(path)
	print('checking ' + path)
	print(str(nAvail) + ' time steps available')
	print(str(nDone)  + ' time steps already done')
	sys.stdout.flush()
	if nAvail > nDone:
		n = nDone
		########################################################################
		# read in files (parallel)
		fileNames   = tools.getFileNames(baseName, n, npc)
		manager     = mp.Manager()
		resultsList = manager.list()
		jobs=[]
		for name in fileNames:
			proc=mp.Process(target=processFile, args=(path,name,resultsList))
			jobs.append(proc)
		nRead=0
		while nRead < npc:
			for n1 in range(myNPC): jobs[nRead+n1].start()
			for n1 in range(myNPC): jobs[nRead+n1].join()
			nRead+=myNPC
			tools.printMem('reading in...')
		tools.printMem('files are read in')
		########################################################################
		# finish up (serial)
		data   = resultsList[0]
		nx     = data.shape[0]
		ny     = data.shape[1]
		nz     = data.shape[2]
		cols   = data.shape[3]
		nTot   = nx*ny*nz*npc
		n1d    = int(np.round(np.power(nTot, 1./3.)))
		res    = np.absolute(data[0,0,0,0]-data[1,0,0,0])
		cLim   = int(res*n1d*100.0); cLim/=200.0;
		cArray = np.arange(-cLim+res/2.0, cLim, res)
		masterArray = np.zeros((n1d, n1d, n1d, cols))
		for data in resultsList:
			xmin = np.amin(data[:,:,:,0])
			ymin = np.amin(data[:,:,:,1])
			zmin = np.amin(data[:,:,:,2])
			xi   = np.argmin(np.absolute(cArray-xmin))
			yi   = np.argmin(np.absolute(cArray-ymin))
			zi   = np.argmin(np.absolute(cArray-zmin))
			masterArray[xi:xi+nx,yi:yi+ny,zi:zi+nz]=data
		########################################################################
		print('writing master file')
		sys.stdout.flush()
		np.save(outDir+baseName+"."+tools.getTimeStepString(n)+".npy", masterArray)
		print('writing is done')
		sys.stdout.flush()
		del masterArray, resultsList
		time.sleep(5)
	else:
		print('all avaliable output is done')
		sys.stdout.flush()
		break
























#
