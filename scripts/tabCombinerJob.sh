#!/bin/bash

#SBATCH -J tabCombiner
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH -o ../../turbTest/run1/bin/out_%j.txt
#SBATCH -e ../../turbTest/run1/bin/error_%j.txt
#SBATCH -p development

# The following commands will be executed when this script is run.
module load python
python tabCombiner.py 64 0 1000 ../../turbTest/run1/bin/ ../../turbTest/run1/bin/combTabs/




#
