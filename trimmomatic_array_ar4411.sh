#!/bin/bash
#SBATCH --array=0-7
#SBATCH --job-name=trimmomatic_ar4411
#SBATCH --output=trim_%A_%a.out
#SBATCH --error=trim_%A_%a.err
#SBATCH -J java2
#SBATCH --cpus-per-task=2
#SBATCH -n 4
#SBATCH --time=06:00:00
#SBATCH --mem=16GB

#this script only performs adaptor trimming function with default settings. Once it is verified that 4 required outputs are being produced, the next script will perform mapping and conversion to bam files.

#loading the required modules in working directory
cd /scratch/ar4411/Homework02
module purge
module load trimmomatic/0.36

# Input: Its the files containing reads, ending in _1.fastq.gz or _2.fastq.gz
FILE1=($(ls *1.fastq.gz))
FILE2=($(ls *2.fastq.gz))

#how to take input files in parallel
INPUT1=${FILE1[$SLURM_ARRAY_TASK_ID]}
INPUT2=${FILE2[$SLURM_ARRAY_TASK_ID]}

#output is taking the whole input, along with the file extension "fastq.gz". Therefore use %%.* to format the files.
OUTPUT1=out_paired_${INPUT1%%.*}.fastq.gz
OUTPUT2=out_unpaired_${INPUT1%%.*}.fastq.gz
OUTPUT3=out_paired_${INPUT2%%.*}.fastq.gz
OUTPUT4=out_unpaired_${INPUT2%%.*}.fastq.gz

#using one set of parameters for all files. More than 80% of the reads are recovered after trimming is done.
#paired files contain more reads / weigh more than unpaired files.

java -jar $TRIMMOMATIC_JAR PE -phred33 $INPUT1 $INPUT2 $OUTPUT1 $OUTPUT2 $OUTPUT3 $OUTPUT4 ILLUMINACLIP:TruSeq3-PE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

#if the script doesn't process ALL the files in parallel due to time or memory constraints, run each file that hasn't been processed, one by one.

#java -jar $TRIMMOMATIC_JAR PE -phred33 SRR1313333_1.fastq.gz SRR1313333_2.fastq.gz out_pairedSRR1313333_1.fastq.gz out_unpairedSRR1313333_1.fastq.gz out_pairedSRR1313333_2.fastq.gz out_unpairedSRR1313333_2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

#java -jar $TRIMMOMATIC_JAR PE -phred33 SRR1313334_1.fastq.gz SRR1313334_2.fastq.gz out_pairedSRR1313334_1.fastq.gz out_unpairedSRR1313334_1.fastq.gz out_pairedSRR1313334_2.fastq.gz out_unpairedSRR1313334_2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

#java -jar $TRIMMOMATIC_JAR PE -phred33 SRR1313335_1.fastq.gz SRR1313335_2.fastq.gz out_pairedSRR1313335_1.fastq.gz out_unpairedSRR1313335_1.fastq.gz out_pairedSRR1313335_2.fastq.gz out_unpairedSRR1313335_2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

#java -jar $TRIMMOMATIC_JAR PE -phred33 SRR1313336_1.fastq.gz SRR1313336_2.fastq.gz out_pairedSRR1313336_1.fastq.gz out_unpairedSRR1313336_1.fastq.gz out_pairedSRR1313336_2.fastq.gz out_unpairedSRR1313336_2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

#java -jar $TRIMMOMATIC_JAR PE -phred33 SRR1313337_1.fastq.gz SRR1313337_2.fastq.gz out_pairedSRR1313337_1.fastq.gz out_unpairedSRR1313337_1.fastq.gz out_pairedSRR1313337_2.fastq.gz out_unpairedSRR1313337_2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

#java -jar $TRIMMOMATIC_JAR PE -phred33 SRR1313338_1.fastq.gz SRR1313338_2.fastq.gz out_pairedSRR1313338_1.fastq.gz out_unpairedSRR1313338_1.fastq.gz out_pairedSRR1313338_2.fastq.gz out_unpairedSRR1313338_2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

#java -jar $TRIMMOMATIC_JAR PE -phred33 SRR1313339_1.fastq.gz SRR1313339_2.fastq.gz out_pairedSRR1313339_1.fastq.gz out_unpairedSRR1313339_1.fastq.gz out_pairedSRR1313339_2.fastq.gz out_unpairedSRR1313339_2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

#java -jar $TRIMMOMATIC_JAR PE -phred33 SRR1313340_1.fastq.gz SRR1313340_2.fastq.gz out_pairedSRR1313340_1.fastq.gz out_unpairedSRR1313340_1.fastq.gz out_pairedSRR1313340_2.fastq.gz out_unpairedSRR1313340_2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

