#!/bin/bash
#SBATCH --job-name=RfBRkrna
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=5-00:00:00
#SBATCH --error=err.braker.rna.log
#SBATCH --output=out.braker.rna.log

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
proteindb=/home/scratch/tpsylvst/databases/Arthropod-protein-db/proteins.fasta
cpu=${SLURM_CPUS_PER_TASK}
software_DIR=/home/scratch/tpsylvst/software
project_DIR=/home/scratch/tpsylvst/projects/$species
genome=$project_DIR/repeatmask-RM/RepeatMasker-cLib/$species.nuclear.purged.fa.masked
script_DIR=$(pwd)
email=tpsylvst@memphis.edu

# rna mappings
file1=$project_DIR/RNAseq-mapping/female/DDM6461/DDM6461.sorted.bam
file2=$project_DIR/RNAseq-mapping/female/DDM6463/DDM6463.sorted.bam
file3=$project_DIR/RNAseq-mapping/male/DDM6462/DDM6462.sorted.bam

# make working dir
mkdir -p $project_DIR/genome-annotation/braker-rna
cd $project_DIR/genome-annotation/braker-rna

# run BRAKER annotation with protein evidence
$software_DIR/BRAKER2/scripts/braker.pl \
--genome=${genome} \
--workingdir=$(pwd) \
--cores $cpu \
--softmasking \
--bam=${file1},${file2},${file3} \
--species=${species}_transcriptome \
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
cat augustus.hints.aa.busco/s* | mail -s "$species braker transcriptome round busco update" $email

# change into the main dir
cd $project_DIR/genome-annotation/braker-protein
