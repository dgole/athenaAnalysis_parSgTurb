#!/bin/bash

#SBATCH -J 3d120
##SBATCH --time=48:00:00
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
##SBATCH -p normal
#SBATCH -p development

# The following commands will be executed when this script is run.
module load gcc python3

export DIR=../../data/turbTest/run30/bin/
python3 quickLook.py $DIR
python3 plots3d.py   $DIR
python3 pspec.py     $DIR
python3 acf.py       $DIR

#export DIR=../../data/turbTest/run30/bin/
#python3 quickLook.py $DIR
#python3 plots3d.py   $DIR
#python3 pspec.py     $DIR
#python3 acf.py       $DIR

#export DIR=../../data/turbTest/run30/bin/
#python3 quickLook.py $DIR
#python3 plots3d.py   $DIR
#python3 pspec.py     $DIR
#python3 acf.py       $DIR

#export DIR=../../data/turbTest/run30/bin/
#python3 quickLook.py $DIR
#python3 plots3d.py   $DIR
#python3 pspec.py     $DIR
#python3 acf.py       $DIR




#
