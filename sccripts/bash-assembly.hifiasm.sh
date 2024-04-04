#!/bin/bash
#SBATCH --job-name=Rfun-asm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=wholeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --error=err.assembly.hifiasm.log
#SBATCH --output=out.assembly.hifiasm.log

# set variables
species=$(ls ../raw-reads/ | sed 's/.fq//')
input=../raw-reads/$species.fq
output=../assembly-hifiasm/$species.hifiasm
cpu=40
scriptDir=$(pwd)

# make output dir
mkdir ../assembly-hifiasm/

# run assembly
/home/tpsylvst/software/hifiasm/hifiasm \
-o $output \
-t $cpu \
$input

# change working dir
cd ../assembly-hifiasm/

# get primary contigs
awk '/^S/{print ">"$2;print $3}' $species.hifiasm.bp.p_ctg.gfa > $species.hifiasm.fa

# get assembly statistics
/home/tpsylvst/software/seqkit/seqkit stats $species.hifiasm.fa > $species.hifiasm.fa.stats

# get busco statistics
# load busco module
module load busco/3.8.7

# return to home dir
cd $scriptDir

# make busco output folder
mkdir ../assembly-hifiasm-busco

# run busco
busco -i ../assembly-hifiasm/$species.hifiasm.fa \
--out_path ../assembly-hifiasm-busco/ \
--out $species \
-m geno \
-c $cpu \
-l endopterygota_odb10
