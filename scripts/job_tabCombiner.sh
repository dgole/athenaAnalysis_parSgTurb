#!/bin/bash

#SBATCH -J tc34
#SBATCH --time=48:00:00
##SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH -p normal
##SBATCH -p development

# The following commands will be executed when this script is run.

export DIR=../../turbTest/run34/bin/

module load gcc python3

cd $DIR
mkdir 1d
mv *.1d ./1d/
cd ../../../athenaAnalysis_parSgTurb/scripts/

echo running quickLook.py
python3 quickLook.py $DIR

#echo running tabCombiner.py
#python3 tabCombiner.py 64 0 1000 $DIR

echo running acf_quick.py
python3 acf_quick.py $DIR 16 6 2



#
