#!/bin/bash
#SBATCH --job-name=mitoHIFI
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G
#SBATCH --time=3-00:00:00
#SBATCH --error=err.mitohifi.log
#SBATCH --output=out.mitohifi.log

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# load conda env
conda activate /home/scratch/tpsylvst/conda-env/mitohifi

# set variables
species=$(ls ../raw-reads/ | sed 's/.fq//')
spDir=$(cd ../ && pwd)
scriptDir=$(pwd)
reads=$scriptDir/../raw-reads/$species.fq
cpu=${SLURM_CPUS_PER_TASK}
MitoHIFI=/home/scratch/tpsylvst/software/MitoHiFi

# set paths
# assembly=/home/scratch/tpsylvst/projects/$species/assembly-hifiasm/$species.hifiasm.fa

# get mitochondria path
Tcas_mito=/home/scratch/tpsylvst/rawdata/fasta/Aromia_bungii_mitochondria/NC_053714.1.fasta
Tcas_mito_gb=/home/scratch/tpsylvst/rawdata/fasta/Aromia_bungii_mitochondria/NC_053714.1.gb

# make output dirs
mkdir -p ../mitochondria/mitohifi
cd ../mitochondria/mitohifi

mkdir reads
cd reads

#minimap2 \
#-t $cpu \
#--secondary=no \
#-ax map-pb \
#/home/scratch/tpsylvst/rawdata/fasta/Aromia_bungii_mitochondria/NC_053714.1.fasta \
#/home/scratch/tpsylvst/projects/Rosalia_funebris/scripts/../raw-reads/Rosalia_funebris.fq | samtools view -@ $cpu -b -F4 -o reads.HiFiMapped.bam

# run mitohifi
python $MitoHIFI/src/mitohifi.py \
-r $reads \
-f $Tcas_mito \
-g $Tcas_mito_gb \
-t $cpu \
-a animal \
--mitos

#cd ..
#mkdir assembly
#cd assembly

# run mitohifi
#python $MitoHIFI/src/mitohifi.py \
#-c $assembly \
#-f $Tcas_mito \
#-g $Tcas_mito_gb \
#-t $cpu \
#-a animal \
#--mitos



