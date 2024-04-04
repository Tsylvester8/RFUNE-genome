#!/bin/bash
#SBATCH --job-name=RfIPR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --error=err.interproscan.log
#SBATCH --output=out.interproscan.log

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# load conda env
conda activate /home/scratch/tpsylvst/conda-env/InterProScan

# this command will remove sequences with internal stop codons and remove ambiguos signs

# set variables
species=$(ls ../raw-reads/ | sed 's/.fq//')
bin_DIR=/home/scratch/tpsylvst/bin
project_DIR=/home/scratch/tpsylvst/projects/$species
script_DIR=$(pwd)
cpu=${SLURM_CPUS_PER_TASK}

# change working dir
mkdir ../genome-annotation/interproscan
cd ../genome-annotation/interproscan
input=$project_DIR/genome-annotation/TSEBRA/TSEBRA.aa

# index fasta file
samtools faidx $input

# get the lines with more than one stop codon
awk -F\* '{print NF-1}' $input > internal-stop-codons.txt

# get protein names
cat $input | grep \> > prot.names.txt

# filter proteins
Rscript $bin_DIR/Rscript.internal.stop.codons.R

# if there are internal stop codons then do the following
if [ -f keep.prots.txt ]; then
    # print file exist
    echo "Removing internal stop codons"
    
    # get final protein list
    samtools faidx $input $(cat keep.prots.txt) > filtered.aa;
    
    # remove terminal stop codons
    cat filtered.aa | sed 's/\*//' > $species.IPR.fasta;
    
    # load busco
    module load busco/3.8.7;
    
    # run busco
    busco -i $species.IPR.fasta -o $species.IPR.busco -m prot -c 10 -l endopterygota_odb10 -f;
         else
         cp $input ./$species.IPR.fasta;
         sed -i 's/\*//' $species.IPR.fasta;
fi

# set interproscan output dir
Interpro_out=$(pwd)

# sort by size
seqkit sort -l $species.IPR.fasta > $species.IPR.sort.fasta

# run InterProScan
/home/scratch/tpsylvst/software/InterProScan/interproscan.sh \
-appl TIGRFAM,FunFam,SFLD,PANTHER,Gene3D,Hamap,PRINTS,Coils,SUPERFAMILY,SMART,CDD,PIRSR,AntiFam,Pfam,MobiDBLite,PIRSF \
-i $species.IPR.sort.fasta \
-f GFF3 \
-cpu $cpu \
-dp \
-goterms

