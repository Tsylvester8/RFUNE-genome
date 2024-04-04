#!/bin/bash
#SBATCH --job-name=Rfmito.map
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --error=err.contaminant.screen.minimap.log
#SBATCH --output=out.contaminant.screen.minimap.log

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# load python conta env
conda activate /home/tpsylvst/conda-env/Genotools/

# load minimap2
module load minimap2/2.17

# set variables
scriptDir=$(pwd)
species=Rosalia_funebris
genome=$scriptDir/../assembly/$species.nuclear.5Mb.fa
reads=../raw-reads/$species.fq
outputDir=$scriptDir/../contaminant-screen
cpu=${SLURM_CPUS_PER_TASK}

# make required dirs
mkdir -p $outputDir
mkdir $outputDir/map

# run minimap2 to map raw reads
minimap2 \
-x map-pb \
-a \
-t $cpu \
-o $outputDir/$species.5Mb.sam \
$genome $reads

# sort mapped reads
samtools sort -@ $cpu -O BAM -o $outputDir/$species.5Mb.sorted.bam $outputDir/$species.5Mb.sam
samtools index $outputDir/$species.5Mb.sorted.bam

# move old files
rm $outputDir/$species.5Mb.sam
