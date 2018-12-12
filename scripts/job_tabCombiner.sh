#!/bin/bash

#SBATCH -J tc12
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH -p normal

# The following commands will be executed when this script is run.

#export DIR=../../turbTest/run34/bin/
export DIR=$SCRATCH/parSgTurb/data/newBCs_turbScaling/run12

module load gcc
module load python3

#cd $DIR
#mkdir 1d
#mv *.1d ./1d/
#cd ../../../athenaAnalysis_parSgTurb/scripts/

#echo running quickLook.py
#python3 quickLook.py $DIR

echo running tabCombiner.py
python3 tabCombiner.py 64 0 1000 $DIR

#echo running acf_quick.py
#python3 acf_quick.py $DIR 16 6 2



#
