#!/bin/bash
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2gb
#SBATCH -o Run_U_5.out
#SBATCH -e Run_U_5.err

../../meiosis -u 1e-05 Run_U_5