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
dvBase   = np.sqrt(alpha)
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/plots3d/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
plt.figure(0)
################################################################################
key = 'dpar'
for n in range(nStart, nEnd):
	data    = do3d.get3d(key, n)
	thisTot = np.sum(data)
	if n = n0:
		tot0 = thisTot
	print(n, str(thisTot/tot0))










#
