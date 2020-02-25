#!/bin/bash
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2gb
#SBATCH -o Run_V_8.out
#SBATCH -e Run_V_8.err

../../meiosis -v 1e-08 Run_V_8