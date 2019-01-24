#!/bin/bash

#SBATCH -J 3d
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 64
#SBATCH -p development

# The following commands will be executed when this script is run.
#module load gcc python3

#export DIR=$SCRATCH/parSgTurb/data/newBCs_turbScaling/run22/bin/
export DIR=../../data/fullPhysicsTest/run64/
python3 plots3d.py   $DIR
python3 pspec.py     $DIR
python3 acf.py       $DIR
python3 slices3d.py  $DIR






#
