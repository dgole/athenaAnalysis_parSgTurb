#!/usr/bin/python
import numpy as np
import matplotlib as m
import scipy
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator
#m.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import math
import sys
sys.path.append('../python')
import athenaReader3d as reader3d
import athenaTools as tools
from matplotlib.backends.backend_pdf import PdfPages
################################################################################
pathBase = str(sys.argv[1])
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/plots3d/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
################################################################################
n = 30
data   = do3d.get3d('dpar', n)
dataCp = data
coordsList = []
coordsList.append((999,999,999))
coordsListRejects = []
for n1 in range(200):
    argLin = np.argmax(dataCp)
    args3d = np.unravel_index(argLin, dataCp.shape)
    print('####################################################')
    print(args3d)
    print(data[args3d])
    dataCp[args3d] = 0.0
    closeFlag = False
    for coords in coordsList:
        dist  = np.square(float(args3d[0]-coords[0]))
        dist += np.square(float(args3d[1]-coords[1]))
        dist += np.square(float(args3d[2]-coords[2]))
        dist  = np.sqrt(dist)
        if dist < 2:
            print("very close to another peak, not recording")
            print(args3d)
            print(coords)
            print(dist)
            closeFlag = True
            coordsListRejects.append(args3d)
            break
    for coords in coordsListRejects:
        dist  = np.square(float(args3d[0]-coords[0]))
        dist += np.square(float(args3d[1]-coords[1]))
        dist += np.square(float(args3d[2]-coords[2]))
        dist  = np.sqrt(dist)
        if dist < 2:
            print("very close to a rejected peak, not recording")
            print(args3d)
            print(coords)
            print(dist)
            closeFlag = True
            coordsListRejects.append(args3d)
            break
    if closeFlag != True: coordsList.append(args3d)
coordsList = coordsList[1:]
print(coordsList)
print(len(coordsList))
#print(coordsListRejects)

massList = []
for coords in coordsList:
    print("assiging mass to coord location")
    gs     = 10
    gsDist = gs*do3d.dx
    i  = coords[0]
    j  = coords[1]
    k  = coords[2]
    localGrid  = data[i-gs:i+gs+1, j-gs:j+gs+1, k-gs:k+gs+1]
    localCoords = np.linspace(0,)
    interpFunc = RegularGridInterpolator(
                (do3d.x[do3d.nx//2-gs:do3d.nx//2+gs+1],
                 do3d.y[do3d.ny//2-gs:do3d.ny//2+gs+1],
                 do3s:k+gs+1]),
                localGrid)
    pts = []
    print(interpFunc(pts))
    print(data[i,j,k])
    print(localGrid[gs,gs,gs])
    mass = np.sum(data[i-2:i+3, j-2:j+3, k-2:k+3])
    if mass > 0.0: massList.append(mass)

nplan   = len(massList)
minMass = min(massList)
sum     = 0
for mass in massList:
    sum += np.log(mass / minMass)
    p   = 1 + nplan * np.power(sum, -1)
    err = (p-1)/np.sqrt(nplan)

print(p,err)







#lp = scipy.ndimage.filters.laplace(data)
#lp /= np.amax(np.absolute(lp))
#plt.imshow(lp[:,:,do3d.nz//2]); plt.colorbar(); plt.show();
#plt.imshow(data[:,:,do3d.nz//2]); plt.colorbar(); plt.show();
#tools.saveAndClear(pathSave + 'timeEvoPertNorm_' + key + '.png', figNum=0)

















#
