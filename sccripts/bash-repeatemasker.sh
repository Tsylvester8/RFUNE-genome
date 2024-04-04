#!/bin/bash
#SBATCH --job-name=RepMask
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=90G
#SBATCH --time=10-00:00:00
#SBATCH --error=err.repeatemasker.log
#SBATCH --output=out.repeatemasker.log

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# activate EDTA environment
conda activate /home/scratch/tpsylvst/conda-env/Repmodeler

# set variables
species=Rosalia_funebris
spDir=$(cd ../ && pwd)
scriptDir=$(pwd)
cpu=${SLURM_CPUS_PER_TASK}
pDIR=/home/scratch/tpsylvst/projects
software=/home/scratch/tpsylvst/software
assembly=$pDIR/$species/assembly/$species.nuclear.purged.fa
ColeoLib=/home/scratch/tpsylvst/databases/Coleoptera-repeat-library/RM-coleoptera-lib.fa
gsize=$(seqkit stats $assembly | awk '{if (NR!=1) print $5}')

## Make working directory
mkdir ../repeatmask-RM
cd ../repeatmask-RM

Rmod_dir=$(pwd)

## buld database for RepeatModeler
$software/RepeatModeler/BuildDatabase \
-name $species \
$assembly

## run RepeatModeler
$software/RepeatModeler/RepeatModeler \
-database $species \
-genomeSampleSizeMax 814001396 \
-threads $cpu \
-LTRStruct > run.out

mkdir RepeatMasker-cLib
cd RepeatMasker-cLib

# make the custom library
cat $ColeoLib ../${species}-families.fa > $species.custom.lib.fa

# run RepeatMasker
$software/RepeatMasker/RepeatMasker \
-e rmblast \
-pa 2 \
-s \
-a \
-xsmall \
-gff \
-dir $(pwd) \
-lib $species.custom.lib.fa \
$assembly

$software/RepeatMasker/util/calcDivergenceFromAlign.pl -s $species.divsum  $species.nuclear.fa.align
tail $species.divsum -n72 > $species.divsum.table
$software/RepeatMasker/util/createRepeatLandscape.pl -div $species.divsum -g 814001396 > $species.html
