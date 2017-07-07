#!/usr/bin/python

####################################################
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import fftpack as fft
import sys
import tabReader as lib
from matplotlib.backends.backend_pdf import PdfPages

dataType=str(sys.argv[1])
runID=int(sys.argv[2])
tsmin=int(sys.argv[3])
tsmax=int(sys.argv[4])
cutFactor=int(sys.argv[5])

outDirName="/home/dan/ResearchLocal/adiabatic/run"+str(runID)+"/"
pp = PdfPages(outDirName+"acf.pdf")

lib.readCombTabs(dataType, runID, tsmin, tsmax, outDirName, cutFactor)


# calculate 2d acf for a single, constant z slice of the data
def acf2dSlice(data, shiftMax, zIndex):
	print "calculating 2d autocorrelation for zIndex = " + str(zIndex)
	corrArray = np.zeros([2*shiftMax+1,2*shiftMax+1])
	# the original 2d slice of data
	array1 = data[...,...,zIndex]
	# loop over all desired x and y shifts
	for i in range(-shiftMax, shiftMax+1):
		#progress = 100*(i+shiftMax)/(2*shiftMax)
		#sys.stdout.write("progress: %d%%\r"%progress)
		#sys.stdout.flush()
		for j in range(-shiftMax, shiftMax+1):
			xShift = i
			yShift = j
			# the 2d slice of data shifted by xShift and yShift
			array2 = np.roll(np.roll(data[...,...,zIndex], xShift, axis=0), yShift, axis=1)
			# loop over all rows in the tables and add the correlations together
			tot = 0
			for k in range(0, lib.numy):
				# correlate one row of the original and shifted arrays
				corr = np.correlate(array1[...,k], array2[...,k])
				tot = tot + corr
			# assign the final sum to it's place in the corrArray based on the x and y shifts
			corrArray[i+shiftMax,j+shiftMax] = tot
	# normalize
	corrArray = corrArray/np.amax(corrArray)
	return corrArray

def getMeanAcf(data, shiftMax, zIndex1, zIndex2):
	for ts in range(0,int((tsmax-tsmin)/cutFactor)):
		for zIndex in range(zIndex1,zIndex2):
			data=lib.bx[ts,...,...,...]
			acf=acf2dSlice(data, shiftMax, zIndex)
			if ts==0 and zIndex==zIndex1:
				acfMean=np.zeros_like(acf)
			acfMean=acfMean+acf
	acfMean=acfMean/float(int((tsmax-tsmin)/cutFactor)+1)/float(zIndex2-zIndex1+1)
	return acfMean
	

acfMeanBx=getMeanAcf(lib.bx, int(lib.numx/2), (lib.numz/2)-5, (lib.numz/2)+5)
acfMeanBy=getMeanAcf(lib.by, int(lib.numx/2), (lib.numz/2)-5, (lib.numz/2)+5)
acfMeanBz=getMeanAcf(lib.bz, int(lib.numx/2), (lib.numz/2)-5, (lib.numz/2)+5)
acfMeanTot=(acfMeanBx+acfMeanBy+acfMeanBz)/3.0
plt.imshow(acfMeanTot, extent=[-1,1,-1,1], cmap=plt.get_cmap("coolwarm"))
plt.xlabel("y/H"); plt.ylabel("x/H");
plt.colorbar(); plt.title("time and space averaged 2d acf"); plt.savefig(pp, format='pdf'); plt.clf();



pp.close()








