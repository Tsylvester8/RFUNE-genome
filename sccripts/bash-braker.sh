#!/bin/bash
#SBATCH --job-name=RfBRAKER
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=10-00:00:00
#SBATCH --error=err.braker.log
#SBATCH --output=out.braker.log

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# activate EDTA environment
conda activate /home/scratch/tpsylvst/conda-env/Braker2

# load augustus module
module load augustus/3.3.1

# load busco
module load busco/3.8.7

# set variables
species=$(ls ../raw-reads/ | sed 's/.fq//')
proteindb=/home/scratch/tpsylvst/databases/Arthropod-protein-db/Arthropoda.fa
cpu=${SLURM_CPUS_PER_TASK}
software_DIR=/home/scratch/tpsylvst/software
project_DIR=/home/scratch/tpsylvst/projects/$species
genome=$project_DIR/repeatmask-RM/RepeatMasker-cLib/$species.nuclear.purged.fa.masked
script_DIR=$(pwd)
email=tpsylvst@memphis.edu

# set output directories
ProtHint_out=$project_DIR/genome-annotation/ProtHint/

# make working dir
mkdir -p $ProtHint_out

# run ProtHint and allign protein sequences to fasta
$software_DIR/ProtHint/bin/prothint.py \
--workdir $ProtHint_out \
--threads $cpu \
$genome \
$proteindb

# run braker on ab initio mode
# make working dir
mkdir -p $project_DIR/genome-annotation/braker-abinitio
cd $project_DIR/genome-annotation/braker-abinitio

## run BRAKER annotation with no evidence
$software_DIR/BRAKER2/scripts/braker.pl \
--genome=${genome} \
--workingdir=$(pwd) \
--cores $cpu \
--softmasking \
--AUGUSTUS_ab_initio \
--esmode \
--species=${species}_abinitio_round \
--rounds 5 \
--AUGUSTUS_BIN_PATH=/public/apps/augustus/3.3.1/bin/ \
--AUGUSTUS_CONFIG_PATH=/home/scratch/tpsylvst/software/Augustus/config \
--AUGUSTUS_SCRIPTS_PATH=/home/scratch/tpsylvst/software/Augustus/scripts \
--GENEMARK_PATH=/home/scratch/tpsylvst/software/GeneMark \
--PROTHINT_PATH=/home/scratch/tpsylvst/software/ProtHint/bin \
--GUSHR_PATH=/home/scratch/tpsylvst/software/GUSHR

# run busco
busco -i augustus.ab_initio.aa -o augustus.ab_initio.aa.busco -m prot -c $cpu -l endopterygota_odb10 -f

## email busco updates
cat augustus.ab_initio.aa.busco/s* | mail -s "$species braker abinitio busco update" $email

## change into the main dir
cd $project_DIR/genome-annotation/braker-abinitio

# run braker with protein evidence mode
# make working dir
mkdir -p $project_DIR/genome-annotation/braker-protein
cd $project_DIR/genome-annotation/braker-protein

# run BRAKER annotation with protein evidence
$software_DIR/BRAKER2/scripts/braker.pl \
--genome=${genome} \
--workingdir=$(pwd) \
--cores $cpu \
--softmasking \
--prot_seq=${proteindb} \
--epmode \
--species=${species}_protein \
--rounds 5 \
--AUGUSTUS_BIN_PATH=/public/apps/augustus/3.3.1/bin/ \
--AUGUSTUS_CONFIG_PATH=/home/scratch/tpsylvst/software/Augustus/config \
--AUGUSTUS_SCRIPTS_PATH=/home/scratch/tpsylvst/software/Augustus/scripts \
--GENEMARK_PATH=/home/scratch/tpsylvst/software/GeneMark \
--PROTHINT_PATH=/home/scratch/tpsylvst/software/ProtHint/bin \
--GUSHR_PATH=/home/scratch/tpsylvst/software/GUSHR

# run busco on annotation with protein homology
busco -i augustus.hints.aa -o augustus.hints.aa.busco -m prot -c $cpu -l endopterygota_odb10 -f

# email busco updates
cat augustus.hints.aa.busco/s* | mail -s "$species braker protein round busco update" $email

# change into the main dir
cd $project_DIR/genome-annotation/braker-protein
