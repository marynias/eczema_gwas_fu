#!/bin/bash

#SBATCH --job-name=garfield
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --time=99:59:00
#SBATCH --mem=60G
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out


/mnt/storage/home/qh18484/scratch/garfield/garfield-v2/garfield