# COMPUTATION / MANIPULATION  
### tabCombiner_big.py  
combines the 3d files into one for each time step in 3d folder  
checks to see if data dump is already done first  
parallelizes up to number of cpus in simulation itself  
args: (npc) (pathToDir1) (pathToDir2) (pathToDir3) ...

# PLOTTING GAS QUANTITIES  
### quickLook.py  
makes plots out of 1d gas averages and particle hist fileList  
args: (path to dir with 1d dir inside)  
saves plots in dir quickLook  

### pspec.py  
power spectrum calculations from 3d files  
args: (path to dir with 3d dir inside)  
saves plots in dir pspec  

### plots3d.py  
plots various profiles and time evolution plots from 3d files  
args: (path to dir with 3d dir inside)  
saves plots in dir plots3d  

### slices.py  
plots slices every 10 dumps for most quantities and every dump for dpar    
args: (path to dir with 3d dir inside)  
saves plots in dir slices  

# PLAN OUTPUT ANALYSIS  
### planPlotter.py  
plots number of particles over time and particle mass histogram  
args: (path to dir with planOutput dir inside)  
saves plots in dir plan  

# CLUMP TRACKING  
### clumpTracker.py   
implements the clump tracking algorithm, which compares particle lists dumped from PLAN to track particles from frame to frame  
args: (path to plan output dir) (start time index) (end time index)  
saves output in ((path to plan output dir)/clumpTracking_(start time index)_ (end time index))  

### ct_fits_mcmc.py  
finds the parmeters of maximum likelihood for various models  
args: (runID)  
saves to (basePath)/plots/clumpTracking_mcmc/  
warning: some of my directory structre and run ID scheme is hard coded into this script currently  

### ct_fits_mcmc_bootstrap_parallel.py  
parallelized bootstrapping of the given fit, dumps parameters to an output file to be read later  
args: (path to plan data) (n_start) (n_end) (number of bootstrap re-samplings) (fit key)  
fit key is something like "bpl", "stpl", etc.  
saves output to (plan path)/plots/clumpTracking_mcmc_bootstrap/  

### ct_fits_mcmc_bootstrap_reader.py  
reads the dumped bootstrap parameters, generates median parameters and error bars, outputs a latex format table  
no args  
note: this is pretty slopily coded right now  

















#
