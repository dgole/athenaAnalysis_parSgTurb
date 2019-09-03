
#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
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
nStart   = int(sys.argv[2])
nEnd     = int(sys.argv[3])
kStart   = int(sys.argv[4])
kEnd     = int(sys.argv[5])
path3d   = pathBase + '3d/'
pathSave = pathBase + '/pspecData/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d     = reader3d.Data3d(path3d)
plt.figure(0)
################################################################################
psk_vy, freqs = reader3d.psProfile_ztAvg(do3d, 'rootRhoDvy', nStart, nEnd, kStart, kEnd)
psk_vx, freqs = reader3d.psProfile_ztAvg(do3d, 'rootRhoDvx', nStart, nEnd, kStart, kEnd)
psk_vz, freqs = reader3d.psProfile_ztAvg(do3d, 'rootRhoDvz', nStart, nEnd, kStart, kEnd)
psk  = psk_vx  + psk_vy  + psk_vz
#psk *= np.power(freqs, -eExpo)
#psk /= np.mean(psk)

saveArr = np.asarray([freqs, psk])
print(saveArr)
fileName = "pspecData_" + str(nStart)+"_"+str(nEnd)+"_"+str(kStart)+"_"+str(kEnd)+".npy"
print("saving file " + fileName)
np.save(pathSave + fileName, saveArr)






#
