# Terrence Sylvester
# pradakshanas@gmail.com
# 25 January 2024

# load library
library(taxize)

# species
sp <- getwd()
sp <- unlist(strsplit(sp,split = "/",fixed = T))[4]
sp <- gsub("_"," ",sp)

# read bacterial and fungal blast hits
PK_blast <- read.delim("../data/hgt/prokaryotic.pblast.out",
                       as.is = T,
                       header = F)
# read animal blast hits
EK_blast <- read.delim("../data/hgt/eukaryotic.pblast.out",
                       as.is = T,
                       header = F)

# read coleoptera blast hits
CL_blast <- read.delim("../data/hgt/coleoptera.pblast.out",
                       as.is = T,
                       header = F)

# rename columns
colnames(PK_blast) <- colnames(EK_blast) <- colnames(CL_blast) <- c("qseqid",
                                                                    "sseqid",
                                                                    "pident",
                                                                    "length",
                                                                    "mismatch",
                                                                    "gapopen",
                                                                    "gaps",
                                                                    "qstart",
                                                                    "qend",
                                                                    "sstart",
                                                                    "send",
                                                                    "evalue",
                                                                    "bitscore",
                                                                    "qcovhsp",
                                                                    "scovhsp")

# get coleoptera taxanomu information
CL_blast$ssp <- NA

for(i in 1:nrow(CL_blast)){
  hit <- unlist(strsplit(CL_blast$sseqid[i], split = "[", fixed = T))
  CL_blast$ssp[i] <- hit[length(hit)]
}
CL_blast$ssp <- gsub("]","",CL_blast$ssp,fixed = T)
CL_blast$ssp <- gsub(":"," ",CL_blast$ssp,fixed = T)

CL_taxonomy <- classification(unique(CL_blast$ssp),db = "ncbi")
CL_blast$ssubfamily <- NA
CL_blast$sfamily <- NA

for(i in 1:length(unique(CL_blast$ssp))){
  loc <- which(names(CL_taxonomy) %in% unique(CL_blast$ssp)[i])
  hit_subfamily <- (CL_taxonomy[[loc]][,1][CL_taxonomy[[loc]][,2] %in% "subfamily"])
  hit_family <- (CL_taxonomy[[loc]][,1][CL_taxonomy[[loc]][,2] %in% "family"])
  
  if(length(hit_subfamily) == 0){
    CL_blast$ssubfamily[CL_blast$ssp %in% unique(CL_blast$ssp)[i]] <- (CL_taxonomy[[loc]][,1][CL_taxonomy[[loc]][,2] %in% "family"])
    CL_blast$sfamily[CL_blast$ssp %in% unique(CL_blast$ssp)[i]] <- (CL_taxonomy[[loc]][,1][CL_taxonomy[[loc]][,2] %in% "family"])
  }else{
    CL_blast$ssubfamily[CL_blast$ssp %in% unique(CL_blast$ssp)[i]] <- hit_subfamily
    CL_blast$sfamily[CL_blast$ssp %in% unique(CL_blast$ssp)[i]] <- hit_family
  }
}

sp_taxonomy <- classification("Rosalia",db = "ncbi")[[1]]

# Neoclytus acuminatus is a Cerambycinae
# Remove all Cerambycinae from coleoptera blast
CL_blast <- CL_blast[CL_blast$ssubfamily != sp_taxonomy$name[sp_taxonomy$rank %in% "subfamily"],]

# read annotation
anno <- read.delim("../data/hgt/Rosalia_funebris.IPR.func.anno.gff3", 
                   as.is = T,
                   header = F,
                   sep = "\t")

# process annotation keep only needed columns
anno <- anno[,c(1,2,9)]
anno <- anno[anno[,2] != ".",]
# rename columns
colnames(anno) <- c("gene",
                    "database",
                    "domain")

for(i in 1:nrow(anno)){
  anno$domain[i] <- unlist(strsplit(unlist(strsplit(anno$domain[i], split = "signature_desc="))[2],
                                    split = ";", fixed = T))[1]
}
# remove empty domains
anno <- anno[!is.na(anno$domain),]

# read braker annotation
Braker_anno <- read.delim("../data/hgt/TSEBRA.gtf", 
                          as.is = T,
                          header = F,
                          sep = "\t",comment.char = "#")

# if Present remove blast hits with an evalue greater than 1e-5
PK_blast <- PK_blast[PK_blast$evalue < 1e-5,]
EK_blast <- EK_blast[EK_blast$evalue < 1e-5,]

# Get blast  hits present only in bacterial and fungal 
PK_only <- PK_blast[!(PK_blast$qseqid %in% EK_blast$qseqid),]

PK_only_tab <-  as.data.frame(matrix(data = NA,
                                     nrow = length(unique(PK_only$qseqid)),
                                     ncol = 10))

colnames(PK_only_tab) <- c("contig",
                           "start",
                           "end",
                           "geneID",
                           "exons",
                           "pID-microbe",
                           "pID-animal",
                           "evalue-microbe",
                           "evalue-animal",
                           "domains")

PK_only_genes <- unique(PK_only$qseqid)

for(i in 1:length(PK_only_genes)){
  tmpMicrobe <- PK_only[PK_only$qseqid %in% PK_only_genes[i],]
  tmpAnno <- Braker_anno[grep(PK_only_genes[i], Braker_anno$V9),]
  PK_only_tab$contig[i] <- unique(tmpAnno$V1)
  PK_only_tab$start[i] <- tmpAnno$V4[tmpAnno$V3 %in% "transcript"]
  PK_only_tab$end[i] <- tmpAnno$V5[tmpAnno$V3 %in% "transcript"]
  PK_only_tab$geneID[i] <- PK_only_genes[i]
  PK_only_tab$exons[i] <- sum(tmpAnno$V3 %in% "exon")
  PK_only_tab$domains[i] <- paste(unique(anno$domain[anno$gene %in% PK_only_genes[i]]), collapse = ", ")
  if(PK_only_tab$domains[i] == ""){
    PK_only_tab$domains[i] <- NA
  }
  
  if(nrow(tmpMicrobe) == 1){
    PK_only_tab$`pID-microbe`[i] <- tmpMicrobe$pident
    PK_only_tab$`evalue-microbe`[i] <- tmpMicrobe$evalue
  }else{
    PK_only_tab$`pID-microbe`[i] <- tmpMicrobe$pident[which.min(tmpMicrobe$evalue)]
    PK_only_tab$`evalue-microbe`[i] <- tmpMicrobe$evalue[which.min(tmpMicrobe$evalue)]
  }
}

# check which PK only hits are present in coleoptera blast
SP_only <- PK_only[!(PK_only$qseqid %in% CL_blast$qseqid),]

SP_only_tab <-  as.data.frame(matrix(data = NA,
                                     nrow = length(unique(SP_only$qseqid)),
                                     ncol = 10))

colnames(SP_only_tab) <- c("contig",
                           "start",
                           "end",
                           "geneID",
                           "exons",
                           "pID-microbe",
                           "pID-animal",
                           "evalue-microbe",
                           "evalue-animal",
                           "domains")

SP_only_genes <- unique(SP_only$qseqid)

for(i in 1:length(SP_only_genes)){
  tmpMicrobe <- PK_only[PK_only$qseqid %in% SP_only_genes[i],]
  tmpAnno <- Braker_anno[grep(SP_only_genes[i], Braker_anno$V9),]
  SP_only_tab$contig[i] <- unique(tmpAnno$V1)
  SP_only_tab$start[i] <- tmpAnno$V4[tmpAnno$V3 %in% "transcript"]
  SP_only_tab$end[i] <- tmpAnno$V5[tmpAnno$V3 %in% "transcript"]
  SP_only_tab$geneID[i] <- SP_only_genes[i]
  SP_only_tab$exons[i] <- sum(tmpAnno$V3 %in% "exon")
  SP_only_tab$domains[i] <- paste(unique(anno$domain[anno$gene %in% SP_only_genes[i]]), collapse = ", ")
  if(SP_only_tab$domains[i] == ""){
    SP_only_tab$domains[i] <- NA
  }
  
  if(nrow(tmpMicrobe) == 1){
    SP_only_tab$`pID-microbe`[i] <- tmpMicrobe$pident
    SP_only_tab$`evalue-microbe`[i] <- tmpMicrobe$evalue
  }else{
    SP_only_tab$`pID-microbe`[i] <- tmpMicrobe$pident[which.min(tmpMicrobe$evalue)]
    SP_only_tab$`evalue-microbe`[i] <- tmpMicrobe$evalue[which.min(tmpMicrobe$evalue)]
  }
}

# Get blast hits present in both bacterial and fungal and animal comparisons
PK_overlap_hits <- PK_blast[(PK_blast$qseqid %in% EK_blast$qseqid),]
EK_overlap_hits <- EK_blast[(EK_blast$qseqid %in% PK_blast$qseqid),]

# get the genes that are more similar to a bacterial or fungal protein than to an
# animal protein

hits <- unique(PK_overlap_hits$qseqid)
putative_HT_loci <- c()
conseved_HT_loci <- c()
for(i in 1:length(hits)){
  tmp_pr <- PK_overlap_hits[PK_overlap_hits$qseqid %in% hits[i],]
  tmp_eu <- EK_overlap_hits[EK_overlap_hits$qseqid %in% hits[i],]
  
  if((nrow(tmp_pr) < 1 | nrow(tmp_eu) < 1) == T){
    stop("We have a problem")
  }
  
  # if the e value of bacterial and fungal blast hit is lower than the
  # e value of animal blast hit we put that in the putative HT list
  if(min(tmp_pr$evalue) < min(tmp_eu$evalue)){
    putative_HT_loci <- c(putative_HT_loci, hits[i])
  }
  
  # if the e value of bacterial and fungal blast hit is higher than the
  # e value of animal blast hit we put that in the conserved HT list
  if(min(tmp_pr$evalue) > min(tmp_eu$evalue)){
    conseved_HT_loci <- c(conseved_HT_loci, hits[i])
  }
}

# get domain names of bacterial and fungal only hits
unique(sort(anno$domain[anno$gene %in% PK_only$qseqid]))

# keep only the gene with highst percentage identity
PK_only_hits <- as.data.frame(matrix(data = NA,
                                     nrow = length(unique(PK_only$qseqid)),
                                     ncol = 4))
colnames(PK_only_hits) <- c("gene",
                            "pident",
                            "nExons",
                            "domain")
PK_only_hits$gene <- unique(PK_only$qseqid)

for(i in 1:nrow(PK_only_hits)){
  tmp <- PK_only[PK_only$qseqid %in% PK_only_hits$gene[i],]
  if(nrow(tmp) == 1){
    PK_only_hits$pident[i] <- PK_only$pident[PK_only$qseqid %in% PK_only_hits$gene[i]]
    PK_only_hits$nExons[i] <- sum(Braker_anno$V3[grep(PK_only_hits$gene[i], Braker_anno$V9)] %in% "exon")
    PK_only_hits$domain[i] <- paste(anno$domain[anno$gene %in% PK_only_hits$gene[i]], collapse = ", ")
  }else{
    PK_only_hits$pident[i] <- max(tmp$pident)
    PK_only_hits$nExons[i] <- sum(Braker_anno$V3[grep(PK_only_hits$gene[i], Braker_anno$V9)] %in% "exon")
    PK_only_hits$domain[i] <- paste(anno$domain[anno$gene %in% PK_only_hits$gene[i]], collapse = ", ")
  }
}

unique(sort(anno$domain[anno$gene %in% conseved_HT_loci]))

putative_HT_loci_tab <-  as.data.frame(matrix(data = NA,
                                              nrow = length(putative_HT_loci),
                                              ncol = 10))

colnames(putative_HT_loci_tab) <- c("contig",
                                    "start",
                                    "end",
                                    "geneID",
                                    "exons",
                                    "pID-microbe",
                                    "pID-animal",
                                    "evalue-microbe",
                                    "evalue-animal",
                                    "domains")

for(i in 1:length(putative_HT_loci)){
  tmpMicrobe <- PK_blast[PK_blast$qseqid %in% putative_HT_loci[i],]
  tmpAnimal <- EK_blast[EK_blast$qseqid %in% putative_HT_loci[i],]
  tmpAnno <- Braker_anno[grep(putative_HT_loci[i], Braker_anno$V9),]
  putative_HT_loci_tab$contig[i] <- unique(tmpAnno$V1)
  putative_HT_loci_tab$start[i] <- tmpAnno$V4[tmpAnno$V3 %in% "transcript"]
  putative_HT_loci_tab$end[i] <- tmpAnno$V5[tmpAnno$V3 %in% "transcript"]
  putative_HT_loci_tab$geneID[i] <- putative_HT_loci[i]
  putative_HT_loci_tab$exons[i] <- sum(tmpAnno$V3 %in% "exon")
  putative_HT_loci_tab$domains[i] <- paste(unique(anno$domain[anno$gene %in% putative_HT_loci[i]]), collapse = ", ")
  if(putative_HT_loci_tab$domains[i] == ""){
    putative_HT_loci_tab$domains[i] <- NA
  }
  
  if(nrow(tmpMicrobe) == 1){
    putative_HT_loci_tab$`pID-microbe`[i] <- tmpMicrobe$pident
    putative_HT_loci_tab$`evalue-microbe`[i] <- tmpMicrobe$evalue
  }else{
    putative_HT_loci_tab$`pID-microbe`[i] <- tmpMicrobe$pident[which.min(tmpMicrobe$evalue)]
    putative_HT_loci_tab$`evalue-microbe`[i] <- tmpMicrobe$evalue[which.min(tmpMicrobe$evalue)]
  }
  if(nrow(tmpAnimal) == 1){
    putative_HT_loci_tab$`pID-animal`[i] <- tmpAnimal$pident
    putative_HT_loci_tab$`evalue-animal`[i] <- tmpAnimal$evalue
  }else{
    putative_HT_loci_tab$`pID-animal`[i] <- tmpAnimal$pident[which.min(tmpAnimal$evalue)]
    putative_HT_loci_tab$`evalue-animal`[i] <- tmpAnimal$evalue[which.min(tmpAnimal$evalue)]
  }
}

write.table(x = PK_only_tab,
            file = "microbe_only_hits.tsv",
            sep = "\t",
            quote = F,
            row.names = F,col.names = T)

write.table(x = putative_HT_loci_tab,
            file = "putative_HT_loci.tsv",
            sep = "\t",
            quote = F,
            row.names = F,col.names = T)

write.table(x = SP_only_tab,
            file = "SP_only.tsv",
            sep = "\t",
            quote = F,
            row.names = F,col.names = T)


