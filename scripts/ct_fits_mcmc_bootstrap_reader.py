#!/usr/bin/python
import numpy as np
import time
#import matplotlib as m
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import athenaTools as tools
import planOutputReader as readerPlan
import ctReader as readerCt
import plan_stats as pstats
################################################################################
# paths and CL args
pathBaseList  = [
				 "../../data/prodRuns/run100/",
				 "../../data/prodRuns/run103/",
				 "../../data/prodRuns/run101/",
				]
runNameList   = ["Control", "Weak", "Moderate"]
nbList        = [500,500,500,500,500,500,100]
dicList = []
list          = [
				[544.9117833528999, 564.2975384531306],
				[59.600078608296144, 69.29295615841158],
				[304.7255560130286, 314.41843356314405],
				[83.35364174733138, 83.35364174733138],
				[0.0, 0.0],
				[128.00649741215818, 128.00649741215818],
				[20.485032812564214, 1.0992777123333326],
				]
dicList.append(list)
list          = [
				[140.0880191100514, 154.58791036318848],
				[47.26408354254062, 54.514029169109165],
				[115.03836372428333, 122.28830935085188],
				[7.835224629241736, 7.835224629241736],
				[0.0, 0.0],
				[39.400434088127724, 39.400434088127724],
				[16.28512675677632, 1.7852355036392282],
				]
dicList.append(list)
list          = [
				[89.14713973362757, 101.8446888132101],
				[9.493223639253273, 15.841998179044538],
				[51.13180017006536, 57.480574709856626],
				[9.817581471494094, 9.817581471494094],
				[0.0, 0.0],
				[26.474126474263784, 26.474126474263784],
				[16.49807128396222, 3.800522204379689]
				]
dicList.append(list)
minMassList = [0.00274, 0.0028, 0.00292 ]

def convert_x_to_m(x, minMass):
	return minMass*np.exp(x)
def convert_name_to_m(xname):
	mname = "$M" + xname[2:]
	return mname
################################################################################
fitInfoList = [
			   pstats.fitInfo_spl,
			   pstats.fitInfo_stpl,
			   pstats.fitInfo_tpl,
			   pstats.fitInfo_bcpl,
			   pstats.fitInfo_bpl,
			   pstats.fitInfo_vtpl,
			   pstats.fitInfo_tspl
			   ]

for n in range(len(pathBaseList)):
	for j in range(len(fitInfoList)):
		nb = nbList[j]
		fitInfo   = fitInfoList[j]
		pathSave  = pathBaseList[n] + 'plots/clumpTracking_mcmc_bootstrap/'
		paramsAll = np.load(pathSave+fitInfo.name+"_"+str(nb)+".npy")
		dic1        = dicList[n][j][0]
		dic2        = dicList[n][j][1]
		masterStr = ""
		if n!=0 and fitInfo.name==fitInfoList[0].name: masterStr+="\hline\n"
		if fitInfo.name==fitInfoList[0].name: masterStr += runNameList[n].ljust(9)
		else      							: masterStr += " ".ljust(9)
		masterStr += "    &    "
		masterStr += fitInfo.name.upper().ljust(5)
		masterStr += "    &    "
		masterStr += str(np.round(dic1,1)).ljust(6)
		masterStr += "    &    "
		masterStr += str(np.round(dic2,1)).ljust(6)
		for i in range(len(fitInfo.paramNames)):
			if fitInfo.paramNames[i][0:2] == "$x":
				#print("converting xs to ms")
				paramsAll[:,i] = convert_x_to_m(paramsAll[:,i], minMassList[n])
				#print(paramsAll[i])
			med  =  np.percentile(paramsAll, 50, axis=0)[i]
			up   =  np.percentile(paramsAll, 84, axis=0)[i] - med
			down = -np.percentile(paramsAll, 16, axis=0)[i] + med
			masterStr += "    &    "
			if fitInfo.paramNames[i][0:2] == "$x":
				name = convert_name_to_m(fitInfo.paramNames[i][:-1])
			else:
				name = fitInfo.paramNames[i][:-1]
			masterStr += name.ljust(10)
			masterStr += " = "
			masterStr += str(np.round(med,4)).ljust(7)
			masterStr += "^{+"
			masterStr += str(np.round(up,4)).ljust(7)
			masterStr += "}_{-"
			masterStr += str(np.round(down,4)).ljust(7)
			masterStr += "}$"
		masterStr += " \\"
		masterStr += "\\"
		print(masterStr)



################################################################################
