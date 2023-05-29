#!/bin/bash
#SBATCH --account=def-yergeaue-ab
#SBATCH --mem-per-cpu=64G
#SBATCH --time=0-12:00
module load gcc/9.3.0 r/4.2.2
Rscript de_d42.R