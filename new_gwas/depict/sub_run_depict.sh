#!/bin/bash

#SBATCH --job-name=depict
#SBATCH --partition=mrcieu2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --time=99:59:00
#SBATCH --mem=60G
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out


cd /mnt/storage/home/qh18484/scratch/depict

./depict.py