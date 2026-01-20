# Generate hetmap and bar plots for CARD resutls
# Ayushi sharam, M.sc Microbiology 2024 sep2024

rm(list=ls(all=TRUE)) #This commnad clears the workspace (will remove all variables)



#Read each file. You need a line per file/genome and to number each 
E1 <- read.delim("E1-Soil.txt", sep='\t')
E2 <- read.delim("E2-River.txt", sep='\t')
E3 <- read.delim("E3- Mus musculas.txt", sep='\t')
E4 <- read.delim("E4-Municipal wastewater.txt", sep='\t')
E5 <- read.delim("E5 food.txt", sep='\t')
C1 <- read.delim("C1-Stool:Rectal sawab.txt", sep='\t')
C2 <- read.delim("C2 Blood sample.txt", sep='\t')
C3 <- read.delim("C3-respiratory.txt", sep='\t')
C4 <- read.delim("C4-CSF.txt", sep='\t')
C5 <- read.delim("C5-pus swab.txt", sep='\t')

#make a presence/absence matrix for AMR genes:
# needs qdapTools package. First time use: install.packages("qdapTools")
library(qdapTools) 
#add all your genomes in the same format:
pam_genes <- mtabulate(list(E1$Best_Hit_ARO, 
                            E2$Best_Hit_ARO,
                            E3$Best_Hit_ARO,
                            E4$Best_Hit_ARO,
                            E5$Best_Hit_ARO,
                            C1$Best_Hit_ARO, 
                            C2$Best_Hit_ARO,
                            C3$Best_Hit_ARO,
                            C4$Best_Hit_ARO,
                            C5$Best_Hit_ARO))
#make all values equal to 1
pam_genes <- (pam_genes>0)*1

#sort genes by most common to least common
pam_genes <- pam_genes[,order(colSums(pam_genes))]
#add row names (you can add the names of your strain or prefered identifier)
rownames(pam_genes) <- c('E1','E2','E3','E4','E5','C1','C2','C3','C4','C5')

#barchart for genes:
barplot(colSums(pam_genes),main = NULL, xlab = "Gene", ylab = "Number of genomes",las=2, cex.names =1)
#barchart for genomes:
barplot(sort(rowSums(pam_genes)),main = NULL, xlab = "Genomes", ylab = "Number of genes",las=2, cex.names =0.5)
#heatmaps:
#oredered with genomes as entered:
heatmap(pam_genes, scale = 'none', col=c('white','gray20'), cexCol =0.5,  Rowv = NA, Colv = NA)
#clustered
heatmap(pam_genes, scale = 'none', col=c('white','gray20'), cexCol =0.5,RowSideColors = c('gray10','gray10','gray10','gray10','gray10','red','red','red','red','red'))

#thetable can be exported to plot in different programs (GraphPad prism):
write.csv(pam_genes,"CARD_gene_table.csv")
#clustered heatmap
h <- heatmap(pam_genes, scale = 'none', col=c('white','gray20'), cexCol =0.5, keep.dendro=TRUE)

#analyse properties of the cluster: 
#get cluster ids for strains and genes
rowClusters <- as.hclust(h$Rowv)
strainClusters <- cutree(rowClusters, k=2)
columnClusters <- as.hclust(h$Colv)
geneClusters <- cutree(columnClusters, k=2)


#Example 1. Look at proportions of mechanisms in strains from each cluster
#make a count matrix for mechansism in each strain:
mechstrain <- mtabulate(list(E1$Resistance.Mechanism, 
                             E2$Resistance.Mechanism,
                             E3$Resistance.Mechanism,
                             E4$Resistance.Mechanism,
                             E5$Resistance.Mechanism,
                             C1$Resistance.Mechanism,
                             C2$Resistance.Mechanism,
                             C3$Resistance.Mechanism,
                             C4$Resistance.Mechanism,
                             C5$Resistance.Mechanism))
rownames(mechstrain) <- c('E1','E2','E3','E4','E5','C1','C2','C3','C4','C5')
sc1 <- colSums(mechstrain[which(strainClusters == 1 ),])
sc1 <- sc1/sum(sc1)*100
sc2 <- colSums(mechstrain[which(strainClusters == 2 ),])
sc2 <- sc2/sum(sc2)*100

#manually add mechanism names, to shoerten. Based on the order in mechstrain:
mechNames = c('efflux','permeability','inactivation','target mod', 'mix1','mix2')
library(viridis) #this is a library to generate colour-blind-friendly pallets
#this baplot compares the proportion ge genes with different mechanisms in the two group of strains
barplot(cbind(sc1,sc2), col = viridis(length(mechNames)), legend.text = mechNames, xlim=c(0, 2 + 2),
        args.legend=list(
          x= 2 +2 ,
          y= 100,
          bty = "n"
        ))

#Example 2. Look at mechanisms of genes in each cluster
#First need to get the mechanism for every gene.
#make a large table with the information for all the genes. Add all your genomes:
allGenes <- rbind(E1,E2,E3,E4,E5,C1,C2,C3,C4,C5)
#get unique gene names
allGenes <- allGenes[!duplicated(allGenes[,'Best_Hit_ARO']),]
gcv <- cbind(cluster =geneClusters,mechanism = rep(0,length(geneClusters))) 

for (i in 1:length(geneClusters)){
  gcv[i,2] <- allGenes[which(allGenes$Best_Hit_ARO == rownames(gcv)[i]),'Resistance.Mechanism']
}

gcv<-as.data.frame(gcv)
temp<- as.data.frame(table(gcv))
ggplot(temp, aes(y=Freq, x=cluster, fill = mechanism)) + 
  geom_bar(position="fill",stat="identity") +
  scale_fill_viridis(discrete = T,labels = mechNames) 
