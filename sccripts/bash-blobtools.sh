#!/bin/bash
#SBATCH --job-name=blobplot
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=10-00:00:00
#SBATCH --error=err.blobtools.log
#SBATCH --output=out.blobtools.log

# set variables
species=$(ls ../raw-reads/ | sed 's/.fq//')
outputDir=../contaminant-screen
cpu=${SLURM_CPUS_PER_TASK}
scriptDir=$(pwd)
pDIR=$scriptDir/../

# set inputs
blast=$scriptDir/../contaminant-screen/mmseqs2/out/Rosalia_funebris.mmseqs.out
map=$pDIR/contaminant-screen/map/$species.5Mb.sorted.bam
genome=$pDIR/assembly/$species.nuclear.5Mb.fa
blobdir=/home/scratch/tpsylvst/software/blobtools

# make required dirs
mkdir -p $outputDir/blobtools

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# load python conta env
conda activate /home/scratch/tpsylvst/conda-env/general-env/

# run blobtools
# create blob database
$blobdir/blobtools create \
-i $genome \
-b $map \
-t $blast \
-o $outputDir/blobtools/$species

# Create a view of a blobDB file
$blobdir/blobtools view \
-i $outputDir/blobtools/$species.blobDB.json \
-o $outputDir/blobtools/$species

# inspect the output
 grep '^##' $outputDir/blobtools/$species.blobDB.table.txt ;  grep -v '^##' $outputDir/blobtools/$species.blobDB.table.txt |  column -t -s $'\t' > $outputDir/analysis/$species.blobtools.table

# Create a blobplot
$blobdir/blobtools plot \
-i $outputDir/blobtools/$species.blobDB.json \
-o $outputDir/blobtools/$species \
--notitle \
--format pdf
