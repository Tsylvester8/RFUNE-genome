#!/bin/bash
#SBATCH --job-name=RfSTAR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpsylvst@memphis.edu
#SBATCH --partition=computeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=10-00:00:00
#SBATCH --error=err.rnamap.star.log
#SBATCH --output=out.rnamap.star.log

# set variables
script_DIR=$(pwd)
species=$(ls $script_DIR/../raw-reads/ | sed 's/.fq//')
cpu=${SLURM_CPUS_PER_TASK}
software_DIR=/home/scratch/tpsylvst/software
project_DIR=/home/scratch/tpsylvst/projects/$species
genome=$project_DIR/assembly/$species.nuclear.purged.fa
fasta=$species.nuclear.purged.fa

# make working dir
mkdir -p $script_DIR/../RNAseq-mapping/genome
cd $script_DIR/../RNAseq-mapping
# copy genome into working directory
cp $genome $script_DIR/../RNAseq-mapping/genome/

# create genome index
$software_DIR/STAR/source/STAR \
--runThreadN $cpu \
--runMode genomeGenerate \
--genomeDir $script_DIR/../RNAseq-mapping/genome \
--genomeFastaFiles $script_DIR/../RNAseq-mapping/genome/$fasta \
--genomeSAindexNbases 13

#echo "####"
#echo "genome index generation complete"
#echo "####"

# align reads
for sex in $(ls -1 $project_DIR/RNAseq/);
do
  for i in $(ls -1 $project_DIR/RNAseq/$sex);
  do
  # make output dir
  mkdir -p $script_DIR/../RNAseq-mapping/$sex/$i;
  # run aligner
  $software_DIR/STAR/source/STAR \
  --runThreadN $cpu \
  --twopassMode Basic \
  --genomeDir $script_DIR/../RNAseq-mapping/genome \
  --readFilesIn $project_DIR/RNAseq/$sex/$i/${i}_1.fq.gz $project_DIR/RNAseq/$sex/$i/${i}_2.fq.gz \
  --outFileNamePrefix $script_DIR/../RNAseq-mapping/$sex/$i/${i}. \
  --readFilesCommand zcat;

  echo "####";
  echo "STAR mapping Diaprepes complete";
  echo "####";

  # sort
  samtools sort \
  --threads $cpu \
  -o $script_DIR/../RNAseq-mapping/$sex/$i/${i}.sorted.bam \
  $script_DIR/../RNAseq-mapping/$sex/$i/${i}.Aligned.out.sam;

# remove old sam file
  rm $script_DIR/../RNAseq-mapping/$sex/$i/${i}.Aligned.out.sam

  # index
  samtools index \
  $script_DIR/../RNAseq-mapping/$sex/$i/${i}.sorted.bam;

  echo "####";
  echo "sorting and indexting complete";
  echo "####";
  done;
done
