#!/bin/bash

#SBATCH -J tc30
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
##SBATCH -o ../../turbTest/run30/bin/out_%j.txt
##SBATCH -e ../../turbTest/run30/bin/error_%j.txt
#SBATCH -p normal

# The following commands will be executed when this script is run.
module load python
cd ../../turbTest/run30/bin/
mkdir 1d
mv *.1d ./1d/
python tabCombiner.py 64 0 1000 ../../turbTest/run30/bin/
python quickLook.py ../../turbTest/run30/bin/
python acf_quick.py ../../turbTest/run30/bin/



#
