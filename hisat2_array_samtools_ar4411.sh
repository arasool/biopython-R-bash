#!/bin/bash
#SBATCH --array=0-7
#SBATCH --job-name=hisat2_paired_ar4411
#SBATCH --output=hisat2_p_%A_%a.out
#SBATCH --error=hisat2_p_%A_%a.err
#SBATCH -J hisat2
#SBATCH --cpus-per-task=2
#SBATCH -n 16
#SBATCH --time=06:00:00
#SBATCH --mem=16GB

#loading the required modules for mapping and conversion to .bam files in working directory.
# will be using paired output from trimmomatic script.
cd /scratch/ar4411/Homework02
module purge
module load hisat2/intel/2.0.5
module load samtools/intel


# capture the current input that will be used for first step, and store it in a variable.
# we are only going to align the paired files since unpaired ones are expected to have much lesser reads,
# so they will not be processed.
FILE1=($(ls out_paired*_1.fastq.gz))
FILE2=($(ls out_paired*_2.fastq.gz))

# this is going to assign the variables to file names.
# all input files should start processing in parallel.
INPUT1=${FILE1[$SLURM_ARRAY_TASK_ID]}
INPUT2=${FILE2[$SLURM_ARRAY_TASK_ID]}

#${f%.*} to get rid of previous file extension.
OUTPUT=${INPUT1%%.*}.sam
OUTPUT2=${OUTPUT%.*}.bam
OUTPUT3=${OUTPUT2%.*}.sorted.bam

# both parts of the pair -1 and -2 will produce one complete sam alignment file.
# the output sam files might be labeled as _1.sam but it will have alignment from both pairs.
# outsplice_paired contains a list of  strand specific coordinates for splice variants.

hisat2 --fr --threads 10 --novel-splicesite-outfile outsplice_paired --rna-strandness FR -x mm10/genome -1 $INPUT1 -2 $INPUT2  -S $OUTPUT


# take the sam file from previous step and convert it to bam
samtools view -o $OUTPUT2 -b $OUTPUT
# sorting the bam files
samtools sort -o $OUTPUT3 $OUTPUT2
# indexing the bam files
samtools index $OUTPUT3

# going to use sorted bam files for featureCounts.
