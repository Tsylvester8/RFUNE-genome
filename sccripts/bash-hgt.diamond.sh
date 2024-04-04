#!/bin/bash
#SBATCH --job-name=RosaliaHGT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=acomputeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=05:00:00
#SBATCH --error=err.hgt.diamond.log
#SBATCH --output=out.hgt.diamond.log

# load module spack
module load spack/0.21.0

# load module diamond
module load  diamond/2.1.7-v5qtpzu

# set variables
species=$((cd .. && pwd) | sed "s+$(cd ../../ && pwd)++" | sed 's+/++')
scriptDir=$(pwd)
cpu=${SLURM_CPUS_PER_TASK}

# make working dir
mkdir -p $scriptDir/../hgt
cd $scriptDir/../hgt

# run diamond on prokaryotic proteins
# diamond blastp \
#	--threads $cpu \
#	--verbose \
#	--db /scratch/tpsylvst/database/diamond/NCBI-prokaryotic-proteins/NCBI-prokaryotic-proteins \
#	--out prokaryotic.pblast.out \
#	--evalue 1e-5 \
#	--ultra-sensitive \
#	--query $scriptDir/../genome-annotation/$species.proteins.fasta \
#	--outfmt 6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovhsp scovhsp

# run diamond on prokaryotic proteins
# diamond blastp \
#	--threads $cpu \
#	--verbose \
#	--db /scratch/tpsylvst/database/diamond/NCBI-eukaryotic-proteins/NCBI-eukaryotic-proteins \
#	--out eukaryotic.pblast.out \
#	--evalue 1e-5 \
#	--ultra-sensitive \
#	--query $scriptDir/../genome-annotation/$species.proteins.fasta \
#	--outfmt 6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovhsp scovhsp

# run diamond on coleoptera proteins
diamond blastp \
        --threads $cpu \
        --verbose \
        --db /scratch/tpsylvst/database/diamond/NCBI-coleoptera-proteins/NCBI-coleoptera-proteins \
        --out coleoptera.pblast.out \
        --evalue 1e-5 \
        --ultra-sensitive \
        --query $scriptDir/../genome-annotation/$species.proteins.fasta \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovhsp scovhsp
