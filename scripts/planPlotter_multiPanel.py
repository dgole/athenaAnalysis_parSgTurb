#!/usr/bin/python
import numpy as np
import matplotlib as m
#m.use('Agg')
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import athenaReader1d as reader1d
import athenaReaderPhst as readerPhst
import athenaTools as tools
import planOutputReader as readerPlan
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import numpy.polynomial.polynomial as poly
################################################################################
# paths and CL args
pathBase = str(sys.argv[1])
nTot     = int(sys.argv[2])
pathPlan = pathBase + 'planOutput/'
pathSave = pathBase + 'plots/planAnim/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
################################################################################
# set up data object and other params
# 1    2    3      4    5    6
# npar mass r_hill xcom ycom zcom
doPlan1 = readerPlan.DataPlan(pathPlan, nStart=300, nTot=nTot, nPar=512*512*512, dt=0.1)
p_mle_master   = []
err_mle_master = []
p_fit_master   = []
alpha_spl_master   = []
alpha_stpl_master   = []
xexp_stpl_master   = []
################################################################################
def makeAnimFrame(doPlan, n):
	# calculate the PL slope with the MLE and fit
	doPlot = True
	try:
		mp, dndmp = readerPlan.getDiffMassHist(doPlan, n)
		p_mle, err_mle        = readerPlan.get_p_mle(doPlan, n)
		p_fit, err_fit, c_fit = readerPlan.get_p_fit(doPlan, n)
		p_mle_master.append(p_mle)
		err_mle_master.append(err_mle)
		p_fit_master.append(p_fit)
	except:
		doPlot = False
		p_mle_master.append(0.0)
		err_mle_master.append(0.0)
		p_fit_master.append(0.0)
		alpha_spl_master.append(0.0)
		alpha_stpl_master.append(0.0)
		xexp_stpl_master.append(0.0)
	if doPlot:
		# Setup
		print('saving anim frame for n = ' + str(n))
		fac = 0.6
		fig = plt.figure(figsize=(19*fac, 22*fac), dpi=120)
		ax = []
		ax.append(plt.subplot2grid((8, 4), (0, 0), rowspan=2, colspan=2)) # 0 row 1 first 3
		ax.append(plt.subplot2grid((8, 4), (0, 2), rowspan=2, colspan=2)) # 1
		ax.append(plt.subplot2grid((8, 4), (2, 0), rowspan=2, colspan=2)) # 2
		ax.append(plt.subplot2grid((8, 4), (2, 2), rowspan=2, colspan=2)) # 3
		ax.append(plt.subplot2grid((8, 4), (4, 0), colspan=4, rowspan=2)) # 3 3 bottom plots
		ax.append(plt.subplot2grid((8, 4), (6, 0), colspan=4)) # 4
		ax.append(plt.subplot2grid((8, 4), (7, 0), colspan=4)) # 5

		# scatter plots
		axNum = 0
		masses = doPlan.peakArrayList[n][:, 2]
		xs     = doPlan.peakArrayList[n][:, 4]
		ys     = doPlan.peakArrayList[n][:, 5]
		zs     = doPlan.peakArrayList[n][:, 6]
		sizes  = [np.power(1.e4*mass, 1./2.) for mass in masses]
		ax[axNum].scatter(xs, ys, s=sizes)
		ax[axNum].set_xlabel(r'$r/h$')
		ax[axNum].set_ylabel(r'y/h')
		ax[axNum].set_ylim(-0.1, 0.1)
		ax[axNum].set_xlim(-0.1, 0.1)
		axNum = 1
		ax[axNum].scatter(xs, zs, s=sizes)
		ax[axNum].set_xlabel(r'$r/h$')
		ax[axNum].set_ylabel(r'$z/h$')
		ax[axNum].set_ylim(-0.1, 0.1)
		ax[axNum].set_xlim(-0.1, 0.1)

		# plot histogram
		axNum = 2
		mp = np.asarray(mp); dndmp = np.asarray(dndmp);
		ax[axNum].loglog(mp, dndmp, 'ko', ms=2)
		# plot MLE model slope
		mp1 = np.logspace(-10.0, 0.0, num=50)
		for i in np.arange(-10,10,0.5):
			preFactor = np.power(10,i)
			model = preFactor*np.power(mp1, -p_mle)
			ax[axNum].loglog(mp1, model, color=(0,0,0,0.2), linestyle='--', linewidth=1)
		# plot average fit to histogram
		model = c_fit * np.power(mp1, -p_fit)
		ax[axNum].loglog(mp1, model, color=(0,0,0,1.0), linestyle='-', linewidth=1)
		# plot labels etc.
		ax[axNum].set_xlabel(r'$M_p$')
		ax[axNum].set_ylabel(r'$dN/dM_p$')
		ax[axNum].set_ylim(1.e1,  1.e8)
		ax[axNum].set_xlim(1e-5, 1.e-1)

		# plot cumulative histogram
		# setup
		axNum  = 3
		mp1, ngtm = readerPlan.getCumMassHist(doPlan, n)
		nm       = mp1.shape[0]
		xCoord = doPlan.tMax+0.4
		yScale = doPlan.nClumpsList[-1]
		yBase  = 12.5
		plt.text(xCoord, (yBase-0.0)*yScale, "###### PARAMS ######")
		# do all fits
		means_spl, errsPlus_spl, errsMinus_spl, maxLike_spl = readerPlan.bootstrap(mp1, readerPlan.fit_spl, 1, nb=1)
		means_stpl, errsPlus_stpl, errsMinus_stpl, maxLike_stpl = readerPlan.bootstrap(mp1, readerPlan.fit_stpl, 2, nb=1)
		means_vtpl, errsPlus_vtpl, errsMinus_vtpl, maxLike_vtpl = readerPlan.bootstrap(mp1, readerPlan.fit_vtpl, 3, nb=1)
		try: means_bcpl, errsPlus_bcpl, errsMinus_bcpl, maxLike_bcpl = readerPlan.bootstrap(mp1, readerPlan.fit_bcpl, 3, nb=1)
		except: pass
		try: means_tpl, errsPlus_tpl, errsMinus_tpl, maxLike_tpl = readerPlan.bootstrap(mp1, readerPlan.fit_tpl, 2, nb=1)
		except: pass
		try: means_bpl, errsPlus_bpl, errsMinus_bpl, maxLike_bpl = readerPlan.bootstrap(mp1, readerPlan.fit_bpl, 3, nb=1)
		except: pass
		# information critera stuff
		printIC = True
		try:
			K = 1
			bic_spl = readerPlan.BIC(K, nm, maxLike_spl)
			aic_spl = readerPlan.AIC(K, nm, maxLike_spl)
			K = 2
			bic_stpl = readerPlan.BIC(K, nm, maxLike_stpl)
			aic_stpl = readerPlan.AIC(K, nm, maxLike_stpl)
			K = 3
			bic_vtpl = readerPlan.BIC(K, nm, maxLike_vtpl)
			aic_vtpl = readerPlan.AIC(K, nm, maxLike_vtpl)
			K = 3
			bic_bcpl = readerPlan.BIC(K, nm, maxLike_bcpl)
			aic_bcpl = readerPlan.AIC(K, nm, maxLike_bcpl)
			K = 2
			bic_tpl = readerPlan.BIC(K, nm, maxLike_tpl)
			aic_tpl = readerPlan.AIC(K, nm, maxLike_tpl)
			K = 3
			bic_bpl = readerPlan.BIC(K, nm, maxLike_bpl)
			aic_bpl = readerPlan.AIC(K, nm, maxLike_bpl)

			bic_min = min(bic_spl, bic_stpl, bic_vtpl, bic_bcpl, bic_tpl, bic_bpl)
			aic_min = min(aic_spl, aic_stpl, aic_vtpl, aic_bcpl, aic_tpl, aic_bpl)

			dbic_spl  = bic_spl  - bic_min
			daic_spl  = aic_spl  - aic_min
			dbic_stpl = bic_stpl - bic_min
			daic_stpl = aic_stpl - aic_min
			dbic_vtpl = bic_vtpl - bic_min
			daic_vtpl = aic_vtpl - aic_min
			dbic_bcpl = bic_bcpl - bic_min
			daic_bcpl = aic_bcpl - aic_min
			dbic_tpl  = bic_tpl  - bic_min
			daic_tpl  = aic_tpl  - aic_min
			dbic_bpl  = bic_bpl  - bic_min
			daic_bpl  = aic_bpl  - aic_min
		except:
			printIC = False

		# SPL
		alpha_spl_master.append(means_spl)
		ngtm_spl = readerPlan.P_spl(mp1, means_spl)
		fac = 1.e5
		ax[axNum].loglog(mp1, fac*ngtm, color=(0,0,0,0.2), marker=".", linewidth=5, markersize=2)
		ax[axNum].loglog(mp1, nm*fac*ngtm_spl, color=(0,0,0,1), label='SPL', linestyle='--')
		plt.text(xCoord, (yBase-0.3)*yScale, "SPL:")
		plt.text(xCoord, (yBase-0.5)*yScale, np.round(means_spl,3))
		if printIC: plt.text(xCoord, (yBase-0.7)*yScale, [np.round(dbic_spl,2), np.round(daic_spl,2)])
		# STPL
		ngtm_stpl = readerPlan.P_stpl(mp1, means_stpl)
		fac = 1.e4
		ax[axNum].loglog(mp1, fac*ngtm, color=(1,0,0,0.2), marker=".", linewidth=5, markersize=2)
		ax[axNum].loglog(mp1, nm*fac*ngtm_stpl, color=(1,0,0,1), label='STPL', linestyle='--')
		plt.text(xCoord, (yBase-1.0)*yScale, "STPL:")
		plt.text(xCoord, (yBase-1.2)*yScale, np.round(means_stpl,3))
		if printIC: plt.text(xCoord, (yBase-1.4)*yScale, [np.round(dbic_stpl,2), np.round(daic_stpl,2)])
		# VTPL
		ngtm_vtpl = readerPlan.P_vtpl(mp1, means_vtpl)
		fac = 1.e3
		ax[axNum].loglog(mp1, fac*ngtm, color=(0,1,0,0.2), marker=".", linewidth=5, markersize=2)
		ax[axNum].loglog(mp1, nm*fac*ngtm_vtpl, color=(0,1,0,1), label='VTPL', linestyle='--')
		plt.text(xCoord, (yBase-1.7)*yScale, "VTPL:")
		plt.text(xCoord, (yBase-1.9)*yScale, np.round(means_vtpl,3))
		if printIC: plt.text(xCoord, (yBase-2.1)*yScale, [np.round(dbic_vtpl,2), np.round(daic_vtpl,2)])
		# BCPL
		try:
			ngtm_bcpl = readerPlan.P_bcpl(mp1, means_bcpl)
			fac = 1.e2
			ax[axNum].loglog(mp1, fac*ngtm, color=(0,0,1,0.2), marker=".", linewidth=5, markersize=2)
			ax[axNum].loglog(mp1, nm*fac*ngtm_bcpl, color=(0,0,1,1.0), label='BCPL', linestyle='--')
			plt.text(xCoord, (yBase-2.4)*yScale, "BCPL:")
			plt.text(xCoord, (yBase-2.6)*yScale, np.round(means_bcpl,3))
			if printIC: plt.text(xCoord, (yBase-2.8)*yScale, [np.round(dbic_bcpl,2), np.round(daic_bcpl,2)])
		except: pass
		# TPL
		try:
			ngtm_tpl = readerPlan.P_tpl(mp1, means_tpl)
			fac =1.e1
			ax[axNum].loglog(mp1, fac*ngtm, color=(1,0,1,0.2), marker=".", linewidth=5, markersize=2)
			ax[axNum].loglog(mp1, nm*fac*ngtm_tpl, color=(1,0,1,1), label='TPL', linestyle='--')
			plt.text(xCoord, (yBase-3.1)*yScale, "TPL:")
			plt.text(xCoord, (yBase-3.3)*yScale, np.round(means_tpl,3))
			if printIC: plt.text(xCoord, (yBase-3.5)*yScale, [np.round(dbic_tpl,2), np.round(daic_tpl,2)])
		except: pass
		# BPL
		try:
			ngtm_bpl = readerPlan.P_bpl(mp1, means_bpl)
			fac = 1.e0
			ax[axNum].loglog(mp1, fac*ngtm, color=(0.3,0.7,1,0.2), marker=".", linewidth=5, markersize=2)
			ax[axNum].loglog(mp1, nm*fac*ngtm_bpl, color=(0.3,0.7,1,1.0), label='BPL', linestyle='--')
			plt.text(xCoord, (yBase-3.8)*yScale, "BPL:")
			plt.text(xCoord, (yBase-4.0)*yScale, np.round(means_bpl,3))
			if printIC: plt.text(xCoord, (yBase-4.2)*yScale, [np.round(dbic_bpl,2), np.round(daic_bpl,2)])
		except: pass
		# other plot stuff
		ax[axNum].set_xlabel(r'$M_p$')
		ax[axNum].set_ylabel(r'$N(>M_p)$')
		ax[axNum].set_ylim(1.e-1, 1.e8)
		ax[axNum].set_xlim(1e-5, 1.e-1)
		ax[axNum].legend(prop={'size':6})


		# plot p values over time
		# p_mle
		axNum = 4
		ax[axNum].plot(doPlan.time[:n+1], p_mle_master[:n+1], 'k', linewidth=1)
		ax[axNum].plot(doPlan.time[n], p_mle_master[n], 'ko', markersize=3)
		pArr = np.asarray(p_mle_master); errArr = np.asarray(err_mle_master)
		ax[axNum].fill_between(doPlan.time[:n+1], pArr-errArr, pArr+errArr, color='gray', alpha=0.3)
		# alpha from spl fit
		ax[axNum].plot(doPlan.time[:n+1], alpha_spl_master[:n+1], 'b', linewidth=1)
		ax[axNum].plot(doPlan.time[n], alpha_spl_master[n], 'bo', markersize=3)
		pArr = np.asarray(alpha_spl_master)
		# alpha and xexp from stpl fit
		#ax[axNum].plot(doPlan.time[:n+1], alpha_stpl_master[:n+1], 'r', linewidth=1)
		#ax[axNum].plot(doPlan.time[n], alpha_stpl_master[n], 'ro', markersize=3)
		#ax[axNum].plot(doPlan.time[:n+1], xexp_stpl_master[:n+1], 'r', linewidth=1)
		#ax[axNum].plot(doPlan.time[n], xexp_stpl_master[n], 'ro', markersize=3)
		#pArr = np.asarray(alpha_spl_master)
		# other stuff
		ax[axNum].axhline(1.4, color=(0,0,0,0.2), linestyle='--')
		ax[axNum].axhline(1.7, color=(0,0,0,0.2), linestyle='--')
		ax[axNum].set_ylim(-0.1,3.0)
		ax[axNum].set_xlim(0.0, doPlan.tMax)
		p_mle_str = r'$p_{mle}=$'  + str(np.round(p_mle,   2)) + r'$\pm$' + str(np.round(err_mle, 2))
		p_fit_str = r'$p_{fit}=$'  + str(np.round(p_fit,   2)) + r'$\pm$' + str(np.round(err_fit, 2))
		ax[axNum].set_title(p_mle_str + "    " + p_fit_str)

		# plot mass frac in planetesimals over time
		axNum = 5
		mFracList = []
		for item in doPlan.peakArrayList:
			try:
				mNow = np.sum(item[:,2])
			except:
				mNow = 0.0
			mFracNow = mNow / doPlan.mParTot
			mFracList.append(mFracNow)
		ax[axNum].plot(doPlan.time[n:], mFracList[n:], 'gray', linewidth=1)
		ax[axNum].plot(doPlan.time[:n+1], mFracList[:n+1], 'k', linewidth=2)
		ax[axNum].plot(doPlan.time[n], mFracList[n], 'ko', markersize=5)
		ax[axNum].set_xlabel(r'$t \Omega$')
		ax[axNum].set_ylabel(r'$M_{plan} / M_{par}$')
		ax[axNum].set_xlim(0.0, doPlan.tMax)

		# plot nClumps over time
		axNum = 6
		ax[axNum].plot(doPlan.time[n:], doPlan.nClumpsList[n:], 'gray', linewidth=1)
		ax[axNum].plot(doPlan.time[:n+1], doPlan.nClumpsList[:n+1], 'k', linewidth=2)
		ax[axNum].plot(doPlan.time[n], doPlan.nClumpsList[n], 'ko', markersize=5)
		ax[axNum].set_ylabel(r'$N_{clumps}$')
		ax[axNum].set_xlabel(r'$t \Omega$')
		ax[axNum].set_xlim(0.0, doPlan.tMax)

		# close and save figure
		#plt.title("n="+str(n))
		plt.tight_layout()
		plt.savefig(pathSave + "anim_" + str(n) + ".png", bbox_inches='tight')
		plt.close('all')

################################################################################

for n in range(0, doPlan1.nTot):
	makeAnimFrame(doPlan1, n)











#
