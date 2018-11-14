#!/bin/bash

#SBATCH -J acfQ30
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
##SBATCH -o ../../parhTest/run12/bin/out_%j.txt
##SBATCH -e ../../parhTest/run12/bin/error_%j.txt
#SBATCH -p normal

# The following commands will be executed when this script is run.
module load python
python acf_quick.py ../../turbTest/run30/bin/





#
