#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4gb
#SBATCH -o runU5.out
#SBATCH -e runU5.err

../../meiosis -u 1e-5 runU5

