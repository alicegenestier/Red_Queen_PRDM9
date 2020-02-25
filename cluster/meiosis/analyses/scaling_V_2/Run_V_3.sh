#!/bin/bash
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2gb
#SBATCH -o Run_V_3.out
#SBATCH -e Run_V_3.err

../../meiosis -v 0.001 Run_V_3