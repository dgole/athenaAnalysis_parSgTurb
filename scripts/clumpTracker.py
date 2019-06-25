#!/usr/bin/python
import numpy as np
#import matplotlib as m
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
#import athenaTools as tools
#import planOutputReader as readerPlan
################################################################################
# paths and CL args
pathBase  = str(sys.argv[1])
nStart    = int(sys.argv[2])
nStop     = int(sys.argv[3])
pathPlan  = pathBase
pathSave  = pathBase + 'clumpTracking_' + str(nStart) + '_' + str(nStop) + '/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
################################################################################

def get_timeStepString(i):
	if i > 999:	 zstring = ""
	elif i > 99: zstring = "0"
	elif i > 9:  zstring = "00"
	elif i > -1: zstring = "000"
	return zstring+str(i)

def get_fileNameList(pathBase, n):
	nStr         = get_timeStepString(n)
	dirName      = pathPlan + 'ParList.' + nStr + '/'
	fileNameList = os.listdir(dirName)
	idList       = []
	for fileName in fileNameList:
		idList.append((n, int(fileName[:-4])))
	return dirName, fileNameList

def get_parIdArr(dirName, fileName):
	#print("reading in file " + dirName + fileName)
	thisIdList = []
	inFile = open(dirName+fileName, 'r')
	for line in inFile:
		split = line.split()
		thisIdList.append(int(split[0]))
	return np.asarray(thisIdList)

class Clump:
	def __init__(self, clumpId, parIdArr):
		#print("initializing a clump data structure...")
		self.idList      = [clumpId]
		self.currentPars = parIdArr
	def update(self, clumpId, parIdArr):
		if parIdArr is not None:
			self.currentPars = parIdArr
			self.idList.append(clumpId)

def compare_arrays(a0, a1):
	n0      = a0.shape[0]
	aBoth   = a0
	aBoth   = np.append(a0, a1)
	aUnique = np.unique(aBoth)
	fracDup0 = (aBoth.shape[0] - aUnique.shape[0]) / a0.shape[0]
	fracDup1 = (aBoth.shape[0] - aUnique.shape[0]) / a1.shape[0]
	return fracDup0, fracDup1

splitterIds = []
splitterData = []
def matchOneClump(clump, parIdArrList, clumpIdList, verbose=False):
	global splitterIds
	if verbose:
		print("###########################################################")
		print("matching current clump " + str(clump.idList[-1]) + " to new frame...")
		sys.stdout.flush()
	matchFracList  = []
	matchFrac1List = []
	nNonZero      = 0
	for parIdArr in parIdArrList:
		matchFrac, matchFrac1 = compare_arrays(clump.currentPars, parIdArr)
		matchFracList.append(matchFrac)
		matchFrac1List.append(matchFrac1)
		if matchFrac > 0.0: nNonZero+=1
	matchFracArr   = np.asarray(matchFracList)
	matchFrac1Arr  = np.asarray(matchFrac1List)
	jArgSort      = (np.argsort(matchFracArr))[::-1]
	if verbose: print("number of non-zero match fractions: " + str(nNonZero))
	if nNonZero > 1:
		i=1
		while matchFracArr[jArgSort[i]] > 0.0:
			if verbose or True:
				print("flagging clump as splitter: " + str(clumpIdList[jArgSort[i]]))
				print("fraction of split clump that was in the old clump: " + str(matchFrac1Arr[jArgSort[i]]) + ", " + str(matchFracArr[jArgSort[i]]))
				sys.stdout.flush()
			splitterIds.append(clumpIdList[jArgSort[i]])
			splitterData.append([matchFrac1Arr[jArgSort[i]], matchFracArr[jArgSort[i]]])
			i+=1
	if nNonZero == 0:
		if verbose:
			print("did not find any matching clumps")
			sys.stdout.flush()
		return None, (None, None), None
	else:
		if verbose:
			print("found" + str(nNonZero) + "clumps with any matching particles")
			print("matching to clump with most overlap: " + str(clumpIdList[jArgSort[0]]))
			print("match fraction was: " + str(matchFracSort[0]))
			sys.stdout.flush()
		return jArgSort[0], clumpIdList[jArgSort[0]], parIdArrList[jArgSort[0]]
################################################################################

clumpObjList = []
print("generating clump objects for starting frame")
sys.stdout.flush()
n = nStart
dirName, fileNameList = get_fileNameList(pathBase, n)
for fileName in fileNameList:
	parIdArr = get_parIdArr(dirName, fileName)
	clumpId  = (n, int(fileName[:-4]))
	clumpObjList.append(Clump(clumpId, parIdArr))

for n in range(nStart+1, nStop):
	print("####################################################################################")
	print("matching to frame " + str(n))
	print("reading in all new clump/particle lists...")
	sys.stdout.flush()
	dirName, fileNameList = get_fileNameList(pathBase, n)
	parIdArrList = []
	clumpIdList  = []
	for fileName in fileNameList:
		parIdArr = get_parIdArr(dirName, fileName)
		clumpId  = (n, int(fileName[:-4]))
		parIdArrList.append(parIdArr)
		clumpIdList.append(clumpId)
	print("looping over all already-existing clumps to find matches...")
	sys.stdout.flush()
	claimedClumpsList = []
	nDone = 0; nTot = len(clumpObjList);
	for clump in clumpObjList:
		if np.mod(nDone, 10) == 0: print("clump " + str(nDone) + " of " + str(nTot))
		sys.stdout.flush()
		nDone+=1
		jMax, clumpIdMatch, parIdArrMatch = matchOneClump(clump, parIdArrList, clumpIdList, verbose=False)
		clump.update(clumpIdMatch, parIdArrMatch)
		claimedClumpsList.append(jMax)
	print("making new clump objects for unclaimed clumps...")
	sys.stdout.flush()
	for j in range(len(fileNameList)):
		#print(claimedClumpsList)
		if j not in claimedClumpsList:
			clumpObjList.append(Clump(clumpIdList[j], parIdArrList[j]))

	print("summary of current state:")
	totClumps        = len(clumpObjList)
	newClumps        = 0
	deadClumps       = 0
	liveClumps       = 0
	startedSplitters = 0
	for clump in clumpObjList:
		if clump.idList[-1][0] != n: deadClumps += 1
		else:                        liveClumps += 1
		if clump.idList[0][0]  == n: newClumps  += 1
		if len(splitterIds)>0:
			#print("clump ID is: " + str(clump.idList[0][1]))
			#print("splitter ID list is: " + str(splitterIdsArr))
			if clump.idList[0] in splitterIds:
				print(clump.idList[0])
				startedSplitters+=1
	print("total clumps:      " + str(totClumps))
	print("live clumps:       " + str(liveClumps))
	print("dead clumps:       " + str(deadClumps))
	print("new clumps:        " + str(newClumps))
	print("all splitters:     " + str(len(splitterIds)))
	print("started splitters: " + str(startedSplitters))
	#print(splitterIds)
	sys.stdout.flush()

print("####################################################################################")
print("final summary:")
for clump in clumpObjList:
	print(clump.idList)
sys.stdout.flush()

print("####################################################################################")
print("saving results...")
sys.stdout.flush()
i = 0
for clump in clumpObjList:
	saveArr = np.asarray(clump.idList)
	fileName = pathSave + str(i) + ".npy"
	np.save(fileName, saveArr)
	i+=1
# save splitterIds
saveArr = np.asarray(splitterIds)
fileName = pathSave + "splitterIds.npy"
np.save(fileName, saveArr)
# save splitterData
saveArr = np.asarray(splitterData)
fileName = pathSave + "splitterData.npy"
np.save(fileName, saveArr)




#
