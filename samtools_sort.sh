#!/bin/bash
#SBATCH --array=0-3
#SBATCH --output=samtools_%A_%a.out
#SBATCH --error=samtools_%A_%a.err
#SBATCH -J samtools
#SBATCH -n 4
#SBATCH --time=01:00:00
#SBATCH --mem=16GB

module purge
module load samtools/intel

# capture the output of a command line and store it in a variable
FILES=($(ls *.bam))

# this is going to assign the variables to file names
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
OUTPUT=${INPUT}.sorted.sam

samtools -o $OUTPUT $INPUT