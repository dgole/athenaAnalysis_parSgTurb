####################################################################
COMPUTATION / MANIPULATION
####################################################################
# tabCombiner.py
combines the 3d files into one for each time step in 3d folder
args: npc tsStart tsEnd pathToDir

# tabCombinerJob.sh
slurm template script to submit tabCombiner.py

# acf_fullCalc.py
calculates 2d acf for every z and time, saves as 4d .npy file
args: pathToDir

acfJob.sh
slurm template script to submit acf_fullCalc.py

####################################################################
SINGLE RUNS
####################################################################
# quickLook.py
makes plots out of 1d gas averages and particle hist fileList
args: (path to dir with 1d dir inside)
saves plots in dir quickLook

# acf_quick.py
makes sparse spacetime averaged acf plots
args: (path to dir with 3d dir inside)
saves plots in dir quickAcf

# acf_fullPlot.py
plots data from acf npy nFile


####################################################################
MULTIPLE RUNS
####################################################################
turbAnalysis.py
