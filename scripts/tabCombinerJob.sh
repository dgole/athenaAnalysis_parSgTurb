#!/bin/bash

#SBATCH -J tc12
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH -o ../../parhTest/run12/bin/out_%j.txt
#SBATCH -e ../../parhTest/run12/bin/error_%j.txt
#SBATCH -p development

# The following commands will be executed when this script is run.
module load python
python tabCombiner.py 64 0 1000 ../../parhTest/run12/bin/ 




#
