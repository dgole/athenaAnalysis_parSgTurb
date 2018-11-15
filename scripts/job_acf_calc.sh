#!/bin/bash

#SBATCH -J acf34
#SBATCH --time=48:00:00
##SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH -p normal
##SBATCH -p development

# The following commands will be executed when this script is run.

module load gcc python3

export DIR=../../turbTest/run34/bin/

python3 acf_calc.py $DIR 4 2 2




#
