#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors
import sys
sys.path.append('../python')
import athenaTools as tools
import planOutputReader as readerPlan
import numpy.polynomial.polynomial as poly
import time
import emcee
import scipy.optimize as op
import corner

################################################################################
nwalkers = 100
nsteps   = 1000
nbrute   = 50

def convert_to_x(masses):
	return np.log(masses/np.amin(masses))

def p_to_P(masses, p):
	#print("from inside p_to_P")
	#print(masses.shape)
	#print(p.shape)
	#dm      = np.zeros_like(p)
	#for i in range(dm.shape[0]-1):
		#dm[i] = masses[i+1] - masses[i]
	#dm[-1]  = 0
	#sumThis = dm * p
	sumThis = p
	P       = np.zeros_like(masses)
	for i in range(masses.shape[0]):
		P[i] = np.sum(sumThis[i:])
	P /= P[0]
	return P
################################################################################

def lnlike(params, masses, pFunc, sign):
	pArr = pFunc(params, masses)
	sum  = np.sum(np.log(pArr))
	#if sum > 0:   sum = -1.e5
	if np.isnan(sum): sum = -1.e10
	return sum * sign

def opt_min(masses, pFunc, params0, bounds=None):
	#'''
	result = op.minimize(lnlike,
						 params0,
						 method = "L-BFGS-B",
			             #bounds=bounds,
						 args=(masses, pFunc, -1.0)
						 )
	#print("op.minimize with L-BFGS-B result:")
	#print(result)
	#'''
	'''
	result = op.minimize(lnlike,
						 params0,
						 method = "Nelder-Mead",
						 options = {'xatol': 1.e-5},
						 args=(masses, pFunc, -1.0)
						 )
	#print("op.minimize with NM result:")
	#print(result)
	'''
	'''
	result = op.minimize(lnlike,
						 params0,
						 method = "BFGS",
						 args=(masses, pFunc, -1.0)
						 )
	print("op.minimize with BFGS result:")
	print(result)
	'''
	return result

def opt_min_brute(masses, pFunc, bounds):
	result = op.brute(lnlike,
					  bounds,
					  Ns = nbrute,
					  args=(masses, pFunc, -1.0),
					  full_output = True
					  )
	return result

def lnlike_withPrior(params, masses, pFunc, priorFunc):
	like  = lnlike(params, masses, pFunc, 1.0)
	prior = priorFunc(params, masses)
	if np.isnan(prior+like):
		return -np.inf
	return prior + like

def mcmc_sample(masses, pFunc, priorFunc, params0, width):
	ndim      = len(params0)
	nwalkers1 = nwalkers
	if ndim > 3: nwalkers1 = 1*nwalkers
	pos = [np.asarray(params0) + np.asarray(width)*np.random.randn(ndim) for i in range(nwalkers1)]
	sampler = emcee.EnsembleSampler(nwalkers1, ndim, lnlike_withPrior,
									args=(masses, pFunc, priorFunc))
	sampler.run_mcmc(pos, nsteps)
	return sampler

def cornerPlot(samples, labelList):
	fig = corner.corner(samples,
						labels=labelList,
						quantiles=[0.16, 0.5, 0.84],
						show_titles=True)
	return fig


class Fit_Pipeline:
	def __init__(self, masses, fitInfo, verbose=True, makePlots=True):
		if verbose:	print("\n\nfitting " + fitInfo.name)
		self.fitInfo = fitInfo
		self.minFound = None
		self.ndim = len(fitInfo.params0)
		if verbose:	print("doing large scale mcmc...")
		self.sampler = mcmc_sample(masses,
								   fitInfo.pFunc,
								   fitInfo.priorFunc,
								   fitInfo.params0,
								   fitInfo.widths)
		self.samples = self.sampler.chain[:, 100:, :].reshape((-1, self.ndim))
		if makePlots: self.cp1 = cornerPlot(self.samples, fitInfo.paramNames)
		self.params_mcmc = np.asarray([np.percentile(self.samples, 50, axis=0),
									   np.percentile(self.samples, 84, axis=0),
									   np.percentile(self.samples, 16, axis=0)])
		self.params_mcmc_like = lnlike(self.params_mcmc[0,:], masses, fitInfo.pFunc, 1)
		if verbose: print("median mcmc params and uncerts and likelihood are:")
		if verbose: print(self.params_mcmc)
		if verbose: print(self.params_mcmc_like)

		if verbose:	print("doing optimization with OPT method...")
		self.res_opt    = opt_min(masses, fitInfo.pFunc, self.params_mcmc[0])
		self.params_opt = self.res_opt["x"]
		if self.res_opt.success == True: #and self.res_bfgs.nit > 2:
			if verbose:	print("success, best params: ")
			if verbose:	print(self.params_opt)
		else:
			if verbose: print("opt optimize call failed")
		self.params_opt_like = lnlike(self.params_opt, masses, fitInfo.pFunc, 1)
		'''
		if verbose:	print("doing optimization with brute force method...")
		bounds = []
		for i in range(len(fitInfo.params0)):
			bounds.append((self.params_mcmc[2,i],self.params_mcmc[1,i]))
			print(bounds)
		self.res_brute    = opt_min_bfgs(masses, fitInfo.pFunc, bounds)
		self.params_brute = self.res_brute["x"]
		if self.res_brute.success == True: #and self.res_bfgs.nit > 2:
			if verbose:	print("success, best params: ")
			if verbose:	print(self.params_brute)
		else:
			if verbose: print("brute optimize call failed")
		self.params_brute_like = lnlike(self.params_brute, masses, fitInfo.pFunc, 1)
		'''
		if verbose:
			print("comparison of likelihoods:")
			print("mcmc params:    " + str(self.params_mcmc_like))
			print("opt params:    " + str(self.params_opt_like))
		if self.params_opt_like > self.params_mcmc_like:
			if verbose: print("success, opt did better than mcmc")
		else:
			print("opt failed to find better params than mcmc")
		if verbose: print("fractional differences in best params:")
		paramDiffs = np.absolute(self.params_mcmc[0] - self.params_opt) / np.absolute(fitInfo.params0)
		if verbose: print(paramDiffs)
		if np.any(paramDiffs>0.2):
			print("CONCERNINGLY LARGE PARAM DIFF")
			print(paramDiffs)
			self.params_opt="failed"
		self.bic = readerPlan.BIC(self.ndim, len(masses), self.params_opt_like)
		self.aic = readerPlan.AIC(self.ndim, len(masses), self.params_opt_like)



################################################################################
class Fit_Info:
	def __init__(self, name, pFunc, priorFunc, params0, widths, bounds, paramNames, color):
		self.name       = name
		self.pFunc      = pFunc
		self.priorFunc  = priorFunc
		self.params0    = params0
		self.widths     = widths
		self.bounds     = bounds
		self.paramNames = paramNames
		self.color1     = (color[0],color[1],color[2],1)
		self.color2     = (color[0],color[1],color[2],0.2)

################################################################################
# SPL FITTING
params0_spl = [0.4]
widths_spl  = [0.1]
bnds_spl = [(-5.0, 5.0)]
def p_spl(params, masses):
	alpha,  = params
	xArr   = convert_to_x(masses)
	return alpha * np.exp(-alpha*xArr)
def prior_spl(params, masses):
	if bnds_spl[0][0] < params[0] < bnds_spl[0][1]:
		return 0.0
	else:
		return -np.inf
fitInfo_spl = Fit_Info("spl", p_spl, prior_spl, params0_spl, widths_spl, bnds_spl,
					   [r"$\alpha$"], (0,0,0))

################################################################################
# STPL FITTING
params0_stpl = [0.2, 3.0]
widths_stpl  = [0.1, 0.5]
bnds_stpl = [(-5.0, 5.0), (0.0, 10.0)]
def p_stpl(params, masses):
	alpha, x_exp = params
	xArr   = convert_to_x(masses)
	top    = alpha + np.exp(-x_exp) * np.exp(xArr)
	bottom = np.exp( alpha*xArr + np.exp(-x_exp)*(np.exp(xArr)-1.0))
	p = top/bottom
	return p
def prior_stpl(params, masses):
	xArr   = convert_to_x(masses)
	if (bnds_stpl[0][0] < params[0] < bnds_stpl[0][1] and
		bnds_stpl[1][0] < params[1] < bnds_stpl[1][1] and
		params[1]<=np.amax(xArr)):
		return 0.0
	else:
		return -np.inf
fitInfo_stpl = Fit_Info("stpl", p_stpl, prior_stpl, params0_stpl, widths_stpl, bnds_stpl,
						[r"$\alpha$", r"$x_{exp}$"], (1,0,0))

################################################################################
# TPL FITTING
params0_tpl = [0.4, 3.0]
widths_tpl  = [0.1, 0.5]
bnds_tpl = [(1.e-5, 5.0), (0.1, 10.0)]
def p_tpl(params, masses):
	a, xt  = params
	xArr   = convert_to_x(masses)
	#print(np.amax(xArr))
	if xt<np.amax(xArr) and np.amax(xArr)<10.0:
		return np.ones_like(masses)*1.e-1000
	top    = a * np.exp(-a*xArr)
	bottom = 1.0 - np.exp(-a*xt)
	res  = top/bottom
	fac  = np.where(xArr>xt, 0, 1)
	res *=fac
	return res
def prior_tpl(params, masses):
	xArr   = convert_to_x(masses)
	if (bnds_tpl[0][0] < params[0] < bnds_tpl[0][1] and
		bnds_tpl[1][0] < params[1] < bnds_tpl[1][1] and
		params[1]>=np.amax(xArr)):
		return 0.0
	else:
		return -np.inf
fitInfo_tpl = Fit_Info("tpl", p_tpl, prior_tpl, params0_tpl, widths_tpl, bnds_tpl,
						[r"$\alpha$", r"$x_{tr}$"], (1,0,1))

################################################################################
# BCPL FITTING
params0_bcpl = [0.1, 1.5, 2.0]
widths_bcpl  = [0.1, 0.1, 0.5]
bnds_bcpl = [(1.e-5, 5.0), (1.e-5, 5.0), (0.1, 10.0)]
def p_bcpl(params, masses):
	a1, a2, xb = params
	xArr   = convert_to_x(masses)
	ltxbr  = np.where(xArr <= xb, 1.0, 0.0)
	gtxbr  = np.where(xArr >  xb, 1.0, 0.0)
	ltxbr *= a1 * np.exp(-a1 * xArr)
	gtxbr *= a2 * np.exp((a2-a1)*xb - a2*xArr)
	return gtxbr + ltxbr
def prior_bcpl(params, masses):
	xArr   = convert_to_x(masses)
	if (bnds_bcpl[0][0] < params[0] < bnds_bcpl[0][1] and
		bnds_bcpl[1][0] < params[1] < bnds_bcpl[1][1] and
		bnds_bcpl[2][0] < params[2] < bnds_bcpl[2][1] and
		params[2]<=np.amax(xArr)):
		return 0.0
	else:
		return -np.inf
fitInfo_bcpl = Fit_Info("bcpl", p_bcpl, prior_bcpl, params0_bcpl, widths_bcpl, bnds_bcpl,
						[r"$\alpha_1$", r"$\alpha_2$", r"$x_{br}$"], (0.3,1,0.8))

################################################################################
# BPL FITTING
params0_bpl = [-1.5, 1.5, 2.0]
widths_bpl  = [0.1,  0.1, 0.5]
bnds_bpl = [(-5.0, 5.0), (-5.0, 5.0), (0.1, 10.0)]
def p_bpl(params, masses):
	a1, a2, xb = params
	c0     = np.power((1./a1)+((1./a2)-(1./a1))*np.power(np.exp(xb),-a1),-1)
	xArr   = convert_to_x(masses)
	ltxbr  = np.where(xArr <= xb, 1.0, 0.0)
	gtxbr  = np.where(xArr >  xb, 1.0, 0.0)
	ltxbr *= c0 * np.exp(-a1 * xArr)
	gtxbr *= c0 * np.exp((a2-a1)*xb - a2*xArr)
	return gtxbr + ltxbr
def prior_bpl(params, masses):
	xArr   = convert_to_x(masses)
	if (bnds_bpl[0][0] < params[0] < bnds_bpl[0][1] and
		bnds_bpl[1][0] < params[1] < bnds_bpl[1][1] and
		bnds_bpl[2][0] < params[2] < bnds_bpl[2][1] and
		params[2]<=np.amax(xArr)):
		return 0.0
	else:
		return -np.inf
fitInfo_bpl = Fit_Info("bpl", p_bpl, prior_bpl, params0_bpl, widths_bpl, bnds_bpl,
					   [r"$\alpha_1$", r"$\alpha_2$", r"$x_{br}$"], (0,0,1))

################################################################################
# VTPL FITTING
params0_vtpl = [0.4, 0.4, 2.0]
widths_vtpl  = [0.1, 0.1, 0.5]
bnds_vtpl = [(-5.0, 5.0), (-5.0, 5.0), (0.1, 10.0)]
def p_vtpl(params, masses):
	alpha, beta, x_exp = params
	xArr   = convert_to_x(masses)
	top    = alpha + beta*np.exp(beta*(x_exp +xArr))
	bottom = np.exp( alpha*xArr + np.exp(beta*x_exp)*(np.exp(beta*xArr)-1.0) )
	return top/bottom
def prior_vtpl(params, masses):
	xArr   = convert_to_x(masses)
	if (bnds_vtpl[0][0] < params[0] < bnds_vtpl[0][1] and
		bnds_vtpl[1][0] < params[1] < bnds_vtpl[1][1] and
		bnds_vtpl[2][0] < params[2] < bnds_vtpl[2][1] and
		params[2]<=np.amax(xArr)):
		return 0.0
	else:
		return -np.inf
fitInfo_vtpl = Fit_Info("vtpl", p_vtpl, prior_vtpl, params0_vtpl, widths_vtpl, bnds_vtpl,
					   [r"$\alpha$", r"$\beta$", r"$x_{exp}$"], (0.2,1,0))

################################################################################
# TSPL FITTING
params0_tspl = [-1.5, 1.0, 2.0, 1.0, 3.0]
widths_tspl  = [0.1, 0.1, 0.1,  0.5, 0.5]
bnds_tspl = [(-5.0, 5.0), (-5.0, 5.0), (1.e-5, 5.0), (0.1, 10.0), (0.1, 10.0)]
def p_tspl(params, masses):
	a1, a2, a3, xb1, xb2 = params
	c1     = np.power(
					  (1./a1) +
					  ((1./a2)-(1./a1))*np.power(np.exp(xb1),-a1) +
					  ((1./a3)-(1./a2))*np.power(np.exp(xb1),a2-a1)*np.power(np.exp(xb2),-a2)
					  ,-1)
	xArr    = convert_to_x(masses)
	ltxbr1  = np.where(xArr <= xb1, 1.0, 0.0)
	middle  = np.where(xb1 < xArr, 1.0, 0.0) * np.where(xb2 > xArr, 1.0, 0.0)
	gtxbr2  = np.where(xArr >  xb2, 1.0, 0.0)
	ltxbr1  *= c1 * np.exp(-a1 * xArr)
	middle  *= c1 * np.exp((a2-a1)*xb1 - a2*xArr)
	gtxbr2  *= c1 * np.exp(-a3*xArr) / np.exp((a1-a2)*xb1 + (a2-a3)*xb2)
	return gtxbr2 + ltxbr1 + middle
def prior_tspl(params, masses):
	xArr   = convert_to_x(masses)
	if (bnds_tspl[0][0] < params[0] < bnds_tspl[0][1] and
		bnds_tspl[1][0] < params[1] < bnds_tspl[1][1] and
		bnds_tspl[2][0] < params[2] < bnds_tspl[2][1] and
		bnds_tspl[3][0] < params[3] < bnds_tspl[3][1] and
		bnds_tspl[4][0] < params[4] < bnds_tspl[4][1] and
		params[3] < params[4] and
		params[3]<=np.amax(xArr) and
		params[4]<=np.amax(xArr)):
		return 0.0
	else:
		return -np.inf
fitInfo_tspl = Fit_Info("tspl", p_tspl, prior_tspl, params0_tspl, widths_tspl, bnds_tspl,
					   [r"$\alpha_1$", r"$\alpha_2$", r"$\alpha_3$", r"$x_{br1}$", r"$x_{br2}$"], (1,0.5,0))




























#
