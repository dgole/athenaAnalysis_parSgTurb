####################################################################
COMPUTATION / MANIPULATION
####################################################################
# tabCombiner.py
combines the 3d files into one for each time step in 3d folder
checks to see if data dump is already done first
parallelizes infineitley
args: <npc> <pathToDir1> <pathToDir2> <pathToDir3> ...

####################################################################
SINGLE RUNS
####################################################################
# quickLook.py
makes plots out of 1d gas averages and particle hist fileList
args: <path to dir with 1d dir inside>
saves plots in dir quickLook

# pspec.py
power spectrum calculations from 3d files
args: <path to dir with 3d dir inside>
saves plots in dir pspec

# pspec_timeEvo.py
power spectrum over time calculations from 3d files
args: <path to dir with 3d dir inside>
saves plots in dir pspec

# acf.py
acf calculations from 3d files
args: <path to dir with 3d dir inside>
saves plots in dir acf

# plots3d.py
plots various profiles and time evolution plots from 3d files
args: <path to dir with 3d dir inside>
saves plots in dir plots3d

# slices3d.py
plots slices every 10 dumps for most quantities and every dump for dpar
args: <path to dir with 3d dir inside>
saves plots in dir slices

####################################################################
MULTIPLE RUNS
####################################################################
# multi_turbAnalysis.py
compares level of turbulence in the simulation to alpha_in for multiple runs
args: none, edit the script itself
saves plots a level above the data directories, edit path in script

# multi_plots3d.py

####################################################################
PLAN OUTPUT
####################################################################
# planPlotter.py
plots number of particles over time and particle mass histogram
args: <path to dir with planOutput dir inside>
saves plots in dir plan
















#
