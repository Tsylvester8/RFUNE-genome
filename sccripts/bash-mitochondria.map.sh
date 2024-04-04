#!/bin/bash
#SBATCH --job-name=mito.map
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=wholeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --error=err.mitochondria.map.log
#SBATCH --output=out.mitochondria.map.log

# source conda module
source /cm/shared/public/apps/miniconda/3.7/etc/profile.d/conda.sh

# load python conta env
conda activate /home/tpsylvst/conda-env/Genotools/

# load minimap2
module load minimap2/2.17

# set variables
species=$(ls ../raw-reads/ | sed 's/.fq//')
reads=../raw-reads/$species.fq
spDir=$(cd ../ && pwd)
scriptDir=$(pwd)
cpu=40

# set output dir
mkdir $spDir/mitochondria/map

# get the potential mitochondrial contig
m_contig=../mitochondria/blast/potential.mitochondrial.contigs.txt

/home/tpsylvst/software/seqkit/seqkit grep -f $m_contig $spDir/mitochondria/blast/$species.100kb.fasta > $spDir/mitochondria/map/$species.m.ctg.fasta

# run minimap2 to map raw reads to potantial mitochondrial sequence
minimap2 \
-x map-pb \
-a \
-t $cpu \
-o $spDir/mitochondria/map/$species.mito.minimap.sam \
$spDir/mitochondria/map/$species.m.ctg.fasta $reads

# sort mapped reads
samtools sort -@ $cpu -O BAM -o $spDir/mitochondria/map/$species.mito.sorted.bam $spDir/mitochondria/map/$species.mito.minimap.sam

# get only mapped reads
samtools view -b -F 4 -o $spDir/mitochondria/map/$species.mito.sorted.q30.bam $spDir/mitochondria/map/$species.mito.sorted.bam
samtools index $spDir/mitochondria/map/$species.mito.sorted.q30.bam

# move old files
mv $spDir/mitochondria/map/$species.mito.sorted.bam  $spDir/mitochondria/map/.old.$species.mito.sorted.bam
mv $spDir/mitochondria/map/$species.mito.minimap.sam $spDir/mitochondria/map/.old.$species.mito.minimap.sam

# convert to fastq
samtools bam2fq $spDir/mitochondria/map/$species.mito.sorted.q30.bam > $spDir/mitochondria/map/$species.mitoreads.fq

# run Unicycler
/home/tpsylvst/software/Unicycler/unicycler-runner.py \
-l $spDir/mitochondria/map/$species.mitoreads.fq \
-o $spDir/mitochondria/map/$species.mitochondria.genome.assembly \
-t $cpu \
--mode conservative
