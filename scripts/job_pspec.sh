#!/bin/bash

#SBATCH -J pspec_bigBox
##SBATCH --time=48:00:00
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
##SBATCH -p normal
#SBATCH -p development

# The following commands will be executed when this script is run.
module load gcc python3

export DIR=../../data/turbTest/run210/bin/
python3 pspec.py     $DIR





#
