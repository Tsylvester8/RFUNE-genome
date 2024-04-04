# Terrence Sylvester
# pradakshanas@gmail.com
# June 5, 2023

# Modified from: https://github.com/oushujun/EDTA/issues/92

# load libraries
library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)

# read kimura distance data
KimuraDistance <- read.csv("../data/repeats/Rosalia_funebris.divsum.table",sep=" ")
#add here the genome size in bp
genomes_size=829378751

# get repeat class
for(i in 1:ncol(KimuraDistance)){
  colnames(KimuraDistance)[i] <- unlist(strsplit(colnames(KimuraDistance)[i],split = ".",fixed = T))[1]
}

# subset data based on the repeat class
DNA <- KimuraDistance[,colnames(KimuraDistance) == "DNA"]
DNA_sum <- rowSums(DNA)

LINE <- KimuraDistance[,colnames(KimuraDistance) == "LINE"]
LINE_sum <- rowSums(LINE)

LTR <- KimuraDistance[,colnames(KimuraDistance) == "LTR"]
LTR_sum <- rowSums(LTR)

SINE <- KimuraDistance[,colnames(KimuraDistance) == "SINE"]
SINE_sum <- rowSums(SINE)

Satellite <- KimuraDistance[,colnames(KimuraDistance) == "Satellite"]
Simple_repeat <- KimuraDistance[,colnames(KimuraDistance) == "Simple_repeat"]
Unkown <- KimuraDistance[,colnames(KimuraDistance) == "Unknown"]

# assign data to a new data fram for plotting
reps <- as.data.frame(matrix(nrow = nrow(KimuraDistance),ncol = 8))
colnames(reps) <- c("Div","DNA","LINE","LTR","SINE","Satellite","SSR", "Unkown")
reps$Div <- KimuraDistance$Div
reps$DNA <- DNA_sum
reps$LINE <- LINE_sum
reps$LTR <- LTR_sum
reps$SINE <- SINE_sum
reps$Satellite <- 0
reps$SSR <- Simple_repeat
reps$Unkown <- Unkown

# set colours
col <-c("#d73027",
        "#fc8d59",
        "#fee090",
        "#e0f3f8",
        "#91bfdb",
        "#4575b4",
        "gray")

# calculate the pecentage of the repeat sequences in the genome
kd_melt = melt(reps,id="Div")
kd_melt$norm = kd_melt$value/genomes_size * 100

# plot
rPlot <- ggplot(kd_melt, aes(fill=variable, y=norm, x=Div)) + 
  geom_bar(position="stack", stat="identity",color="black", width= 0.75, linewidth = 0.05) +
  scale_fill_viridis(discrete = T) +
  scale_fill_manual(values=col) +
  theme_classic() +
  xlab("Kimura substitution level (CpG adjusted)") +
  ylab("Percent of the genome") + 
  labs(fill = "") +
  ylim(c(0,6)) +
  coord_cartesian(xlim = c(0, 50)) +
  theme(axis.text=element_text(size=11),
        axis.title =element_text(size=12),
        panel.grid.major.y = element_line(colour = "darkgray",linewidth = 0.05),
        panel.grid.minor.y = element_line(colour = "lightgray",linewidth = 0.05))

rPlot

# save plot
ggsave(
  "repeatLandscape.pdf",
  plot = rPlot,
  device = "pdf",
  path = getwd(),
  scale = 1,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

# plot without unknowns
# assign data to a new data fram for plotting
reps <- as.data.frame(matrix(nrow = nrow(KimuraDistance),ncol = 7))
colnames(reps) <- c("Div","DNA","LINE","LTR","SINE","Satellite","SSR")
reps$Div <- KimuraDistance$Div
reps$DNA <- DNA_sum
reps$LINE <- LINE_sum
reps$LTR <- LTR_sum
reps$SINE <- SINE_sum
reps$Satellite <- 0
reps$SSR <- Simple_repeat

# set colours
col <-c("#d73027",
        "#fc8d59",
        "#fee090",
        "#e0f3f8",
        "#91bfdb",
        "#4575b4",
        "gray")

# calculate the pecentage of the repeat sequences in the genome
kd_melt = melt(reps,id="Div")
kd_melt$norm = kd_melt$value/genomes_size * 100



# plot
rPlot <- ggplot(kd_melt, aes(fill=variable, y=norm, x=Div)) + 
  geom_bar(position="stack", stat="identity",color="black", width= 0.75, linewidth = 0.05) +
  scale_fill_viridis(discrete = T) +
  scale_fill_manual(values=col) +
  theme_classic() +
  xlab("Kimura substitution level (CpG adjusted)") +
  ylab("Percent of the genome") + 
  labs(fill = "") +
  ylim(c(0,6)) +
  coord_cartesian(xlim = c(0, 50)) +
  theme(axis.text=element_text(size=11),
        axis.title =element_text(size=12),
        panel.grid.major.y = element_line(colour = "darkgray",linewidth = 0.05),
        panel.grid.minor.y = element_line(colour = "lightgray",linewidth = 0.05))

rPlot

# save plot
ggsave(
  "repeatLandscape-unkown.pdf",
  plot = rPlot,
  device = "pdf",
  path = getwd(),
  scale = 1,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)
