#!/bin/bash

#SBATCH -J bs_all_100
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH --tasks-per-node 1
#SBATCH -p development

# The following commands will be executed when this script is run.
module load gcc
module load python3

python3 ct_fits_mcmc_bootstrap.py ./planNew2/run_plan_100/ 220 275 100
python3 ct_fits_mcmc_bootstrap.py ./planNew2/run_plan_101/ 300 390 100
python3 ct_fits_mcmc_bootstrap.py ./planNew2/run_plan_103/ 220 290 100








#
