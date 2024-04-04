#!/bin/bash
#SBATCH --job-name=Brand-QM
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --time=3-00:00:00
#SBATCH --error=err.mitochondria.blast.log
#SBATCH --output=out.mitochondria.blast.log

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# set variables
species=$(ls ../raw-reads/ | sed 's/.fq//')
reads=../raw-reads/$species.fq
spDir=$(cd ../ && pwd)
scriptDir=$(pwd)
cpu=40

# make output dirs
mkdir -p ../mitochondria/blast

# set paths
assembly=/home/scratch/tpsylvst/projects/$species/assembly-hifiasm/$species.hifiasm.fa
# get blast db path
Tcas_mitochondria_blastdb=/home/scratch/tpsylvst/databases/BLASTDB/Tcas.Mito/Tcas-mitochondria
Pcer_mitochondria_blastdb=/home/scratch/tpsylvst/databases/BLASTDB/Pcer.mito/Pcer-mitochondria
# get mitochondria path
Pcer_mito=/home/scratch/tpsylvst/rawdata/fasta/Polydrusus/Pcer.mitochondria.fasta
Tcas_mito=/home/scratch/tpsylvst/rawdata/fasta/Tribolium/Tcas-mitochondria.fasta

# set blastdb
mitodb=$Tcas_mitochondria_blastdb
#set mitochondria
mito=$Tcas_mito

# filter fasta file by sequence length
/home/tpsylvst/software/seqkit/seqkit seq -M 100000 $assembly > ../mitochondria/blast/$species.100kb.fasta

# load blast module
module load blast/2.12.0

# run blast
blastn -query ../mitochondria/blast/$species.100kb.fasta \
-db $mitodb \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 10 \
-max_hsps 1 \
-evalue 1e-25 \
-out ../mitochondria/blast/$species.mitochondria.blast.out \
-num_threads $cpu

# get potential mitochondrial contigs
cat ../mitochondria/blast/$species.mitochondria.blast.out | awk '{print $1}' > ../mitochondria/blast/potential.mitochondrial.contigs.txt

/home/tpsylvst/software/seqkit/seqkit grep \
-f ../mitochondria/blast/potential.mitochondrial.contigs.txt \
../mitochondria/blast/$species.100kb.fasta > ../mitochondria/blast/potential.mitochondrial.contigs.fasta

# combine with tcas mitochondria for LAST analysis
cat $mito ../mitochondria/blast/potential.mitochondrial.contigs.fasta > ../mitochondria/blast/$species.MAFFT.LAST.fasta
