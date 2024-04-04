#!/bin/bash
#SBATCH --job-name=RFPDups
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=65G
#SBATCH --time=1-00:00:00
#SBATCH --error=err.purgedups.hifiasm.log
#SBATCH --output=out.purgedups.hifiasm.log

# load python and minimap
module load python/3.9.13/gcc.8.2.0
module load minimap2/2.17

# set vars
scriptDir=$(pwd)
species=Rosalia_funebris
cpu=${SLURM_CPUS_PER_TASK}
assembly=$scriptDir/../assembly/$species.nuclear.fa
reads=$scriptDir/../raw-reads/$species.fq
binDir=/home/tpsylvst/software/purge_dups/bin/

# make dirs
mkdir ../assembly-hifiasm-pdups/

# change working dir
cd ../assembly-hifiasm-pdups/

# Run minimap2 to align pacbio data and generate paf files, then calculate read depth
minimap2 \
-xasm20 \
-t $cpu \
$assembly \
$reads | gzip -c - > $species.paf.gz

$binDir/pbcstat $species.paf.gz
$binDir/calcuts PB.stat > cutoffs 2>calcults.log

# Split an assembly and do a self-self alignment
$binDir/split_fa \
$assembly > $species.hifiasm.fa.split

minimap2 \
-xasm5 \
-DP \
-t $cpu \
$species.hifiasm.fa.split \
$species.hifiasm.fa.split | gzip -c - > $species.hifiasm.fa.split.self.paf.gz

# Purge haplotigs and overlaps
$binDir/purge_dups \
-2 \
-T cutoffs \
-c PB.base.cov \
$species.hifiasm.fa.split.self.paf.gz > dups.bed 2> purge_dups.log

# Get purged primary and haplotig sequences from draft assembly
$binDir/get_seqs \
-e dups.bed $assembly

# rename purged assembly
mv purged.fa $species.nuclear.purged.fa

# get assembly stats
/home/tpsylvst/software/seqkit/seqkit stats -a $species.nuclear.purged.fa > $species.nuclear.purged.fa.stats

# run busco
# load busco module
module load busco/3.8.7

# change directory
cd $scriptDir
mkdir ../assembly-hifiasm-pdups-busco/
cd ../assembly-hifiasm-pdups-busco/

# run busco
busco -i ../assembly-hifiasm-pdups/$species.nuclear.purged.fa \
--out $species.dup.purged \
-m geno \
-c $cpu \
-f \
-l endopterygota_odb10
