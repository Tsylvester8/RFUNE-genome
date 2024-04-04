#!/bin/bash
#SBATCH --job-name=RFmmseqs2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --time=15-00:00:00
#SBATCH --error=err.contaminant.screen.mmseqs2.log
#SBATCH --output=out.contaminant.screen.mmseqs2.log

# set variables
scriptDir=$(pwd)
species=Rosalia_funebris
genome=$scriptDir/../assembly/$species.nuclear.5Mb.fa
outputDir=$scriptDir/../contaminant-screen
cpu=${SLURM_CPUS_PER_TASK}
mmseqs2db=/home/scratch/tpsylvst/databases/BLASTDB/NCBI-NT.mmseqs2db/nt.fnaDB
tmp=$outputDir/mmseqs2/tmp

# make required dirs
mkdir -p $outputDir
mkdir $outputDir/mmseqs2

# make other subdirs
mkdir $outputDir/mmseqs2/tmp
mkdir $outputDir/mmseqs2/queryDB
mkdir $outputDir/mmseqs2/search
mkdir $outputDir/mmseqs2/out

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# activate EDTA environment
conda activate /home/scratch/tpsylvst/conda-env/mmseqs2

# change working directory
cd $outputDir/mmseqs2/queryDB

# make fasta file into a database
mmseqs createdb \
--dbtype 2 \
$genome \
$outputDir/mmseqs2/queryDB/$species.queryDB

# return
cd $scriptDir

# change working dir
cd $outputDir/mmseqs2/search

# run mmseqs search
mmseqs search \
--search-type 3 \
--threads $cpu \
-e 1.000E-25 \
--max-seqs 10 \
$outputDir/mmseqs2/queryDB/$species.queryDB \
$mmseqs2db \
$outputDir/mmseqs2/search/$species.resultDB \
$outputDir/mmseqs2/tmp

# return
cd $scriptDir

# change working dir
cd $outputDir/mmseqs2/out

# convert the output to a table
mmseqs convertalis \
--threads $cpu \
--format-mode 0 \
--format-output query,taxid,bits \
$outputDir/mmseqs2/queryDB/$species.queryDB \
$mmseqs2db \
$outputDir/mmseqs2/search/$species.resultDB \
$outputDir/mmseqs2/out/$species.mmseqs.out
