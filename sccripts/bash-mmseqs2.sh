#!/bin/bash
#SBATCH --job-name=mmseqs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=170G
#SBATCH --time=28-00:00:00
#SBATCH --error=err.mmseqs2.log
#SBATCH --output=out.mmseqs2.log

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# activate EDTA environment
conda activate /home/scratch/tpsylvst/conda-env/mmseqs2

# set variables
species=$(ls ../raw-reads/ | sed 's/.fq//')
spDir=$(cd ../ && pwd)
scriptDir=$(pwd)
cpu=${SLURM_CPUS_PER_TASK}
pDIR=/home/scratch/tpsylvst/projects
assembly=$pDIR/$species/assembly/$species.nuclear.fa
db=/home/scratch/tpsylvst/databases/BLASTDB/NCBI-NT.mmseqs2db/nt.fnaDB
tmp=$pDIR/$species/contaminant-screen/mmseq2/tmp

# make temperory directory
mkdir -p $tmp

# Make working directory
mkdir -p $pDIR/$species/contaminant-screen/mmseq2
cd $pDIR/$species/contaminant-screen/mmseq2

# make fasta file into a database
mkdir -p queryDB
cd queryDB

mmseqs createdb \
--dbtype 2 \
$assembly \
./$species.queryDB

# return
cd $scriptDir

# change working dir
mkdir -p $pDIR/$species/contaminant-screen/mmseq2/search
cd $pDIR/$species/contaminant-screen/mmseq2/search

# run mmseqs search
mmseqs search \
--search-type 3 \
--threads $cpu \
-e 1.000E-25 \
--max-seqs 10 \
$pDIR/$species/contaminant-screen/mmseq2/queryDB/$species.queryDB \
$db \
./$species.resultDB \
$pDIR/$species/contaminant-screen/mmseq2/tmp

# return
cd $scriptDir

# change working dir
mkdir -p $pDIR/$species/contaminant-screen/mmseq2/out
cd $pDIR/$species/contaminant-screen/mmseq2/out

# convert the output to a table
mmseqs convertalis \
--threads $cpu \
--format-mode 0 \
--format-output query,taxid,bits \
$pDIR/$species/contaminant-screen/mmseq2/queryDB/$species.queryDB \
$db \
$pDIR/$species/contaminant-screen/mmseq2/search/$species.resultDB \
./$species.mmseqs.out
