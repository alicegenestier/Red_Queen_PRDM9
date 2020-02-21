#!/bin/bash
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2gb
#SBATCH -o Run_N_6.out
#SBATCH -e Run_N_6.err

../../meiosis -N 1000000.0 Run_N_6