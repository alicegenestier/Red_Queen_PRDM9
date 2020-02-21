#!/bin/bash
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2gb
#SBATCH -o Run_simul.out
#SBATCH -e Run_simul.err

../../meiosis Run_simul