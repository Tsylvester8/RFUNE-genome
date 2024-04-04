#!/bin/bash
#SBATCH --job-name=bam2fq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --error=err.bam2fq.log
#SBATCH --output=out.bam2fq.log

## DDM6463 is Rosalia funebris

# set vars
input=/home/scratch/tpsylvst/rawdata/novogene/usftp21.novogene.com/HiFi_data/DDM6463/DDM6463_m64291e_230228_231127.hifi_reads.bam

# make output dir if not already done
mkdir ../raw-reads

# change working dir
cd ../raw-reads

# convert bam to fastq
samtools bam2fq $input > Rosalia_funebris.fq
