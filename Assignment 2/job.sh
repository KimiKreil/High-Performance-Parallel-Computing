#!/usr/bin/env bash
#SBATCH --job-name=Kimi
#SBATCH --partition=modi_short
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --exclusive

mpiexec singularity exec \
   ~/modi_images/hpc-notebook_latest.sif \
   ./task_farm_HEP