#!/bin/bash
#SBATCH --job-name=AGRFblast
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=bigmemq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=5G
#SBATCH --time=06:00:00
#SBATCH --error=err.blast.comps.log
#SBATCH --output=out.blast.comps.log

# load blast module
module load blast/2.12.0

# set vars
script_Dir=$(pwd)
projects_Dir=/home/scratch/tpsylvst/projects
cpu=${SLURM_CPUS_PER_TASK}

# change dir
cd /home/scratch/tpsylvst/projects/Rosalia_funebris/blast-comparisons

# run blast
blastp \
-num_threads $cpu \
-db /home/scratch/tpsylvst/projects/Rosalia_funebris/blast-comparisons/blastdb/AGRF \
-query AGRF.fa \
-outfmt "6 qseqid sseqid evalue" \
-seg yes > blast_output.abc

# -outfmt 7 \
# -seg yes > blast_output.txt

bash process.blast.out.sh

# make output
# mkdir -p $script_Dir/../blast-comparisons
# cd $script_Dir/../blast-comparisons

# do a self blast to calculate copy numbers for a given gene
# blastp -db $projects_Dir/Rosalia_funebris/blastdb/RF_prot \
# -query $projects_Dir/Rosalia_funebris/genome-annotation/TSEBRA/Longest.isomer.aa \
# -out $script_Dir/../blast-comparisons/Rosalia.aa.selfblast.out \
# -outfmt 6 \
# -max_hsps 1 \
# -num_threads $cpu

# blastp -db $projects_Dir/Anoplophora_glabripennis/blastdb/ALB_prot \
# -query $projects_Dir/Rosalia_funebris/genome-annotation/TSEBRA/Longest.isomer.aa \
# -out $script_Dir/../blast-comparisons/Rosalia.aa.vs.ALB.aa.out \
# -outfmt 6 \
# -max_hsps 1 \
# -num_threads $cpu

# blastp -db $projects_Dir/Anoplophora_glabripennis/blastdb/ALB_prot \
# -query $projects_Dir/Anoplophora_glabripennis/genome-annotation/ALB_proteins.faa \
# -out $script_Dir/../blast-comparisons/ALB.aa.selfblast.out \
# -outfmt 6 \
# -max_hsps 1 \
# -num_threads $cpu
