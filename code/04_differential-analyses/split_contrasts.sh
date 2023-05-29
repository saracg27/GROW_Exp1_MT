#!/bin/bash
#SBATCH --account=def-yergeaue-ab
#SBATCH --mem-per-cpu=128G
#SBATCH --time=0-24:00
module load gcc/9.3.0 r/4.2.2
Rscript Split_contrasts.R
