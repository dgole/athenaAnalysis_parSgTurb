#!/bin/bash

#SBATCH -J 3d_100
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH --tasks-per-node 1
#SBATCH -p development

# The following commands will be executed when this script is run.
module load gcc
module load python3


#DIR=$SCRATCH/parSgTurb/data/prodRuns/run101_noPar/bin/
#python3 plots3d.py $DIR 5 9 3.2e-4

#DIR=$SCRATCH/parSgTurb/data/prodRuns/run102_noPar/bin/
#python3 plots3d.py $DIR 5 9 1.e-3

#DIR=$SCRATCH/parSgTurb/data/prodRuns/run103_noPar/bin/
#python3 plots3d.py $DIR 5 9 1.e-4

DIR=$SCRATCH/parSgTurb/data/prodRuns/run100/bin/
python3 plots3d.py $DIR 16 20 1.e-5 

#DIR=$SCRATCH/parSgTurb/data/prodRuns/run101/bin/
#python3 plots3d.py $DIR 16 20 3.2e-4

#DIR=$SCRATCH/parSgTurb/data/prodRuns/run102/bin/
#python3 plots3d.py $DIR 16 20 1.e-3

#DIR=$SCRATCH/parSgTurb/data/prodRuns/run103/bin/
#python3 plots3d.py $DIR 16 20 1.e-4






#python3 multi_pspec.py $DIR

#python3 master_multiPanel.py 1 1.e-3 0.1 $DIR

#python3 plots3d_par.py 64 $DIR
#python3 plots3d.py   $DIR
#python3 pspec.py     $DIR
#python3 pspec_timeEvo.py $DIR
#python3 acf.py       $DIR
#python3 slices3d.py  $DIR
#python3 gas_multiPanelAnimation.py $DIR 1.e-4







#
