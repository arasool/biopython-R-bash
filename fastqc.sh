#!/bin/sh

#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH --mem=4GB
#SBATCH --nodes=1

cd /scratch/ar4411/AG2017

module load fastqc

fastqc KCL_1.fastq

