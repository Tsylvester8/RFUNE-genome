#!/bin/bash
#SBATCH --job-name=Rfune-kmer
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=04:00:00
#SBATCH --error=err.kmer.analysis.log
#SBATCH --output=out.kmer.analysis.log

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

scriptDir=$(pwd)

for i in {21,25,27,29,31,33};

do

# set vars
CPU=${SLURM_CPUS_PER_TASK}
MEM=50
input=../raw-reads/Rosalia_funebris.fq
kmer=$i
outDir=../genome-size-est/${kmer}-mer/
tmp=../genome-size-est/tmp
# make required directories if not already done
mkdir ../genome-size-est
mkdir $outDir
mkdir $tmp

# run kmer analysis
/home/tpsylvst/software/KMC/bin/kmc -k$kmer -m$MEM -t$CPU -ci1 -cs100000 $input $outDir/reads $tmp

# change working dir
cd $outDir

# process kmers
/home/tpsylvst/software/KMC/bin/kmc_tools transform reads histogram reads.histo -cx100000

# activate conda for genomescope
conda activate /home/scratch/tpsylvst/conda-env/GnomeScope

# run genomescope
/home/scratch/tpsylvst/software/genomescope2.0/genomescope.R -i reads.histo -o $(pwd) -k $kmer

# return
cd $scriptDir;

done
