#!/bin/bash
#SBATCH --array=0-3
#SBATCH --output=Bowtie_%A_%a.out
#SBATCH --error=Bowtie_%A_%a.err
#SBATCH -J Bowtie
#SBATCH -n 4
#SBATCH --time=01:00:00
#SBATCH --mem=16GB

module purge
module load bowtie2/intel

# capture the output of a command line and store it in a variable
FILES=($(ls *.fastq))

# this is going to assign the variables to file names
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
OUTPUT=${INPUT}.sam

bowtie2 --phred64 -x Athaliana -p 4 -U $INPUT -S $OUTPUT


