cafeComp <- read.delim("../data/cafe/unfiltered_cafe_input.txt", header = F)
cafeComp <- cafeComp[,c(-1,-5)]
colnames(cafeComp) <- cafeComp[1,]
cafeComp <- cafeComp[-1,]
rownames(cafeComp) <- cafeComp$`Family ID`
cafeComp$XP <- as.numeric(cafeComp$XP)
cafeComp$Rosalia <- as.numeric(cafeComp$Rosalia)

Rosalia <- cafeComp[,c(1,3)]
ALB <- cafeComp[,c(1,2)]

colnames(Rosalia) <- colnames(ALB) <- c("Orthogroup","Count")
Rosalia$Species <- "Rosalia"
ALB$Species <- "ALB"

cafe <- rbind(Rosalia,ALB)
cafe$Orthogroup <- as.numeric(cafe$Orthogroup)
cafe$Count <- as.numeric(cafe$Count)


library(ggplot2)

# Grouped
ggplot(cafe, aes(fill=Species, y=Count, x=Orthogroup)) + 
  geom_bar(position="dodge", stat="identity",width = 0.85) + 
  scale_fill_discrete(type = c("#66c2a5","#fc8d62")) +
  theme_bw()

RosaliaDiff <- (table(sort(cafeComp$Rosalia[cafeComp$Rosalia > cafeComp$XP & cafeComp$XP != 0] - cafeComp$XP[cafeComp$Rosalia > cafeComp$XP & cafeComp$XP != 0])))
ALBDiff <- table(sort(abs(cafeComp$Rosalia[cafeComp$Rosalia < cafeComp$XP & cafeComp$Rosalia != 0] - cafeComp$XP[cafeComp$Rosalia < cafeComp$XP & cafeComp$Rosalia != 0])))

par(mfcol = c(1,2))

barplot(RosaliaDiff,
        xlab = "Orthogroup copy number difference",
        ylab = "Count",
        ylim = c(0,1000),
        col = "lightblue",las = 2,
        space = 0.1)
mtext(text = "A",side = 3,adj = 0,line = 1,cex = 2)

barplot(ALBDiff,
        xlab = "Orthogroup copy number difference",
        ylab = "Count",
        ylim = c(0,1000),
        col = "lightblue",las = 2,
        space = 0.1)
mtext(text = "B",side = 3,adj = 0,line = 1,cex = 2)


