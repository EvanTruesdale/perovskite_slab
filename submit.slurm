#!/bin/bash

#SBATCH --job-name=pv_slab
#SBATCH --account=cmt
#SBATCH --partition=compute

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=8:00:00
#SBATCH --mem=185G

#SBATCH --mail-type=ALL
#SBATCH --mail-user=evan@truecrew.net

#SBATCH --chdir=.
#SBATCH --export=all

module load intel/oneAPI/2021.1.1
source /mmfs1/home/etrues/miniconda3/etc/profile.d/conda.sh
conda activate /mmfs1/home/etrues/miniconda3

## Usual calculations
python3 section_3_A.py 160 4

