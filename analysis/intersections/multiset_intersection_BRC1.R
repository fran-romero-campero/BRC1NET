# R script for pre-processing 
# Copyright (C) 2018  Francisco J. Romero-Campero, Pedro de los Reyes Rodríguez
# Ana Belén Romero Losada
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public
# License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Authors: Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
#          Francisco J. Romero-Campero
# 
# Contact: Francisco J. Romero-Campero - fran@us.es

# Date: March 2019

# Analysis of intersections among multiple sets is fundamental 
# for in-depth understanding of their complex relationships

# This script performs intersections between the targets of two transcription factors and 
# a set of genes. 

#install.packages("SuperExactTest")
library(SuperExactTest)

# ##Reading the two sets TF target genes
# tf1 <- read.table(file = "../../../web_apps/peak_visualizer/data/targets_in_network/CCA1_ZT02_targets_in_network.txt",
#                            header = FALSE, as.is = TRUE)[[1]]
# tf2 <- read.table(file="../../../web_apps/peak_visualizer/data/targets_in_network/LHY_targets_in_network.txt",
#                           header = FALSE, as.is = TRUE)[[1]]
# 
# #Reading the group of genes peaking at specific time
# genes.peak.zt <- read.table(file = "../../../network/clusters/peak_ZT0.txt",
#                              header = FALSE, as.is = TRUE)[[1]]
# 
# #Establishing the list of sets to test
# sets <- list(tf1, tf2, genes.peak.zt)
# 
# 
# #####Test and visualization of intersections#####
# 
# #vignette("set_html")
# results <- supertest(x = sets, n = 5778)
# help("plot.msets")
# par(mar=c(3,3,3,3))
# plot(results, sort.by = "size")
# plot(results, Layout = "landscape")
# 
# results.table <- summary(results)
# 
# #-- Exploring the output--#
# 
# 
# typeof(results.table)
# length(results.table)
# names(results.table)
# results.table$Barcode
# results.table$otab
# typeof(results.table$otab)
# final.intersection <- results.table$otab[["111"]]
# results.table$etab
# results.table$P.value
# tail(results.table$P.value, n=1)
# enrichment <- (results.table$Table)[["FE"]][nrow(results.table$Table)]
# 
# intersection.genes <- (results.table$Table)[["Elements"]][nrow(results.table$Table)]
# 
# intersection.genes <- strsplit(intersection.genes, split = ", ")[[1]]
# 
# length(sets)
# 
# length.gene.sets <- sapply(X = sets,FUN = length)
# 
# #If you want to carry out the test introducing only numbers:
# cpsets(x = 43 -1, L = length.gene.sets, n = 6830, lower.tail = FALSE)




###### ---- Function that returns p-value, enrichment and the set of genes in the intersecion (intersectSets)#####

##Get translation between AGI and primary symbol
library(org.At.tair.db)
columns(org.At.tair.db)
my.key <- keys(org.At.tair.db, keytype="ENTREZID")
my.col <- c("SYMBOL", "TAIR")
alias2symbol.table <- select(org.At.tair.db, keys=my.key, columns=my.col, keytype="ENTREZID")
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias[is.na(alias)] <- "" 


## Get description of each AGI symbol
network.data <- read.table(file="../../data/brc1_transcriptional_network.tsv", 
                           as.is = TRUE, header = TRUE, sep="\t", quote = "", comment.char = "")
head(network.data)

description <- network.data$S.annotation
names(description) <- network.data$names
description[1:3]

#This is the function, called intersectSets
intersectSets <- function(tf1,tf2,set.of.genes, alias,gene.descriptions){
  intersection.data <- list()
  sets <- list(tf1, tf2, set.of.genes)
  #names(sets) <- c("cca1", "lhy", "peakZT0")
  results <- supertest(x = sets, n = 931)
  results.table <- summary(results)
  p.value <- tail(results.table$P.value, n=1) #Get the last p-value
  enrichment <- (results.table$Table)[["FE"]][nrow(results.table$Table)]
  intersection.genes <- (results.table$Table)[["Elements"]][nrow(results.table$Table)]
  intersection.genes <- strsplit(intersection.genes, split = ", ")[[1]]

  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment

  
  intersection.genes.agi <- intersection.genes
  intersection.genes.primary.symbol <- alias[intersection.genes]
  names(intersection.genes.primary.symbol) <- NULL
  gene.table <- matrix(nrow=length(intersection.genes), ncol=3)
  colnames(gene.table) <- c("AGI", "Primary Symbol", "Description")
  gene.table[,1] <- intersection.genes.agi
  gene.table[,2] <- intersection.genes.primary.symbol
  #  gene.table[,3] <- description
  
  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment
  intersection.data[[3]] <- gene.table
  

  intersection.genes.description <- gene.descriptions[intersection.genes]
  names(intersection.genes.description) <- NULL
  

  intersection.data[[3]] <- data.frame(intersection.genes,intersection.genes.primary.symbol,intersection.genes.description,stringsAsFactors = F)
  
  names(intersection.data) <- c("p-value", "enrichment", "gene.table")
  return(intersection.data)
  
}

intersectSets(tf1 = tf1, tf2 = tf2, set.of.genes = genes.peak.zt, alias=alias,gene.descriptions = description)



### Intersection between nodes (genes) with high topological values and clusters of BRC1 network####

brc1.network.data <- read.table(file="../../data/brc1_transcriptional_network.tsv", 
                               sep = "\t", as.is = TRUE, header = TRUE, quote = "", comment.char = "")
head(brc1.network.data)
nrow(brc1.network.data)
gene.names <- brc1.network.data$names

threshold <- 0.75 #Here you can change the threshold

indegree.threshold <- quantile(brc1.network.data$indegree, prob=threshold)
indegree.top <- gene.names[brc1.network.data$indegree > indegree.threshold]

outdegree.threshold <- quantile(brc1.network.data$outdegree, prob=threshold)
outdegree.top <- gene.names[brc1.network.data$outdegree > outdegree.threshold]

brc1.network.degree <- brc1.network.data$indegree + brc1.network.data$outdegree
degree.threshold <- quantile(brc1.network.degree, prob=threshold)
degree.top <- gene.names[brc1.network.degree > degree.threshold]

brc1.network.data$transitivity[is.na(brc1.network.data$transitivity)] <- 0
trans.threshold <- quantile(brc1.network.data$transitivity, prob=threshold)
trans.top <- gene.names[brc1.network.data$trans > trans.threshold]

closeness.threshold <- quantile(brc1.network.data$closeness, prob=threshold)
closeness.top <- gene.names[brc1.network.data$closeness > closeness.threshold]

betweeness.threshold <- quantile(brc1.network.data$betweeness, prob=threshold)
betweeness.top <- gene.names[brc1.network.data$betweeness > betweeness.threshold]

eccentricity.threshold <- quantile(brc1.network.data$eccentricity, prob=threshold)
eccentricity.top <- gene.names[brc1.network.data$eccentricity > eccentricity.threshold]


##--Function to perform an intersection of TWO sets--##
intersect2sets <- function(set1, set2, alias, gene.descriptions){
  intersection.data <- list()
  sets <- list(set1, set2)
  results <- supertest(x = sets, n = 931)
  results.table <- summary(results)
  p.value <- tail(results.table$P.value, n=1) #Get the last p-value
  enrichment <- (results.table$Table)[["FE"]][nrow(results.table$Table)]
  intersection.genes <- (results.table$Table)[["Elements"]][nrow(results.table$Table)]
  intersection.genes <- strsplit(intersection.genes, split = ", ")[[1]]
  
  intersection.genes.agi <- intersection.genes
  intersection.genes.primary.symbol <- alias[intersection.genes]
  names(intersection.genes.primary.symbol) <- NULL
  gene.table <- matrix(nrow=length(intersection.genes), ncol=3)
  colnames(gene.table) <- c("AGI", "Primary Symbol", "Description")
  gene.table[,1] <- intersection.genes.agi
  gene.table[,2] <- intersection.genes.primary.symbol
  #  gene.table[,3] <- description
  


  
  intersection.genes.description <- gene.descriptions[intersection.genes]
  names(intersection.genes.description) <- NULL
  
  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment
  intersection.data[[3]] <- data.frame(intersection.genes,intersection.genes.primary.symbol,intersection.genes.description,stringsAsFactors = F)
  
  names(intersection.data) <- c("p-value", "enrichment", "gene.table")
  return(intersection.data)
  
}

intersect2sets(set1=degree.top, set2 = genes.peak.zt, alias=alias, gene.descriptions = description)

#####---Several if loops to selectize (in the app) the topological parameter and the set of genes---#####

input <- list(cluster=1, topological_parameter="Degree", threshold="0.90")

if (topological_parameter == "Degree")
{
  brc1.network.degree <- brc1.network.data$indegree + brc1.network.data$outdegree
  degree.threshold <- quantile(brc1.network.degree, prob=input$threshold)
  top.genes <- gene.names[brc1.network.degree > degree.threshold]
} else if (topological_parameter == "Transitivity")
{
  brc1.network.trans <- brc1.network.data$transitivity
  trans.threshold <- quantile(brc1.network.trans, prob= input$threshold)
  top.genes <- gene.names[brc1.network.trans > trans.threshold]
} else if (topological_parameter == "Closeness")
{
  brc1.network.closeness <- brc1.network.data$closeness
  closeness.threshold <- quantile(brc1.network.closeness, prob= input$threshold)
  top.genes <- gene.names[brc1.network.closeness > closeness.threshold]
} else if (topological_parameter == "Betweeness")
{
  brc1.network.bet <- brc1.network.data$betweeness
  bet.threshold <- quantile(brc1.network.bet, prob= input$threshold)
  top.genes <- gene.names[brc1.network.bet > bet.threshold]
} else if (topological_parameter == "Eccentricity")
{
  brc1.network.eccen <- brc1.network.data$eccentricity
  eccen.threshold <- quantile(brc1.network.eccen, prob= input$threshold)
  top.genes <- gene.names[brc1.network.eccen > eccen.threshold]
} 

#Classify a gene in its cluster
strsplit(x = brc1.network.data$T.cluster, split = "|")[[1]][1] #cluster1
strsplit(x = brc1.network.data$T.cluster, split = "|")[[8]][6] #cluster2
strsplit(x = brc1.network.data$T.cluster, split = "|")[[17]][11] #cluster3
strsplit(x = brc1.network.data$T.cluster, split = "|")[[15]][16] #cluster4
strsplit(x = brc1.network.data$T.cluster, split = "|")[[22]][21] #cluster5
strsplit(x = brc1.network.data$T.cluster, split = "|")[[20]][26] #cluster6
strsplit(x = brc1.network.data$T.cluster, split = "|")[[82]][31] #cluster7
strsplit(x = brc1.network.data$T.cluster, split = "|")[[45]][36] #cluster8
strsplit(x = brc1.network.data$T.cluster, split = "|")[[81]][41] #cluster9

#For loop to classify all genes in clusters
brc1.clusters <- rep( list(list()), 9 ) 

for (i in 1:nrow(brc1.network.data))
{
  print(i)
  if (strsplit(x = brc1.network.data$T.cluster[i], split = "|")[[1]][1] == 1)
  {
    brc1.clusters[[1]] <- append(brc1.clusters[[1]], gene.names[i])
  }
  
  if (strsplit(x = brc1.network.data$T.cluster[i], split = "|")[[1]][6] == 2)
  {
    brc1.clusters[[2]] <- append(brc1.clusters[[2]], gene.names[i])
  }
  
  if (strsplit(x = brc1.network.data$T.cluster[i], split = "|")[[1]][11] == 3)
  {
    brc1.clusters[[3]] <- append(brc1.clusters[[3]], gene.names[i])
  }
  
  if (strsplit(x = brc1.network.data$T.cluster[i], split = "|")[[1]][16] == 4)
  {
    brc1.clusters[[4]] <- append(brc1.clusters[[4]], gene.names[i])
  }
  
  if (strsplit(x = brc1.network.data$T.cluster[i], split = "|")[[1]][21] == 5)
  {
    brc1.clusters[[5]] <- append(brc1.clusters[[5]], gene.names[i])
  }
  
  if (strsplit(x = brc1.network.data$T.cluster[i], split = "|")[[1]][26] == 6)
  {
    brc1.clusters[[6]] <- append(brc1.clusters[[6]], gene.names[i])
  }
  
  if (strsplit(x = brc1.network.data$T.cluster[i], split = "|")[[1]][31] == 7)
  {
    brc1.clusters[[7]] <- append(brc1.clusters[[7]], gene.names[i])
  }
  
  if (strsplit(x = brc1.network.data$T.cluster[i], split = "|")[[1]][36] == 8)
  {
    brc1.clusters[[8]] <- append(brc1.clusters[[8]], gene.names[i])
  }
  
  if (strsplit(x = brc1.network.data$T.cluster[i], split = "|")[[1]][41] == 9)
  {
    brc1.clusters[[9]] <- append(brc1.clusters[[9]], gene.names[i])
  }
}

length(brc1.clusters[[1]])
##Selectize cluster (for the app)
if (input$cluster == 1)
{
  cluster <- brc1.clusters[[1]]
} else if (input$cluster == 2)
{
  cluster <- brc1.clusters[[2]]
} else if (input$cluster == 3)
{
  cluster <- brc1.clusters[[3]]  
} else if (input$cluster == 4)
{
  cluster <- brc1.clusters[[4]]  
} else if (input$cluster == 5)
{
  cluster <- brc1.clusters[[5]]  
} else if (input$cluster == 6)
{
  cluster <- brc1.clusters[[6]]  
} else if (input$cluster == 7)
{
  cluster <- brc1.clusters[[7]]  
} else if (input$cluster == 8)
{
  cluster <- brc1.clusters[[8]]  
} else if (input$cluster == 9)
{
  cluster <- brc1.clusters[[9]]  
} 

# intersect2sets(set1 = degree.top, brc1.clusters[[8]], alias = alias, gene.descriptions = description)

result <- intersect2sets(set1 = top.genes, set2 = zt.genes, alias = alias, gene.descriptions = description)
p.value <- result[1][[1]]
enrichment <- result[2][[1]]
intersect.genes <- result[3][[1]]$intersection.genes



##Loop to intersect the clusters of genes and high topological values genes####
top.parameters <- c("Degree","Betweeness", "Closeness", "Eccentricity","Transitivity")
top.genes <- list(degree.top, trans.top, closeness.top, betweeness.top, eccentricity.top)
names(top.genes) <- c("Degree", "Transitivity", "Closeness", "Betweeness", "Eccentricity")
#Initialize matrix to store the results
intersection.table <- matrix(ncol=5, nrow = length(brc1.clusters))
colnames(intersection.table) <- c("cluster", "p-value", "fdr", "enrichment", "Intersection Genes") 
head(intersection.table)
i <- 1
j <- 2    
for (i in 1:length(top.parameters))
{
  for (j in 1:length(brc1.clusters))
  {
    
      current.top <- top.genes[i][[1]]
      set.of.genes <- brc1.clusters[[j]]
      
      
      print("TEST")
      result <- intersect2sets(set1 = current.top, set2 = set.of.genes, alias = alias, gene.descriptions = description)
      p.value <- result[1][[1]]
      enrichment <- result[2][[1]]
      intersect.genes <- result[3][[1]]$intersection.genes
        
      
          
      intersection.table[j,1]<- j
      intersection.table[j,2] <- p.value
      intersection.table[j,4] <- enrichment
      intersection.table[j,5] <- paste(intersect.genes, collapse= ",")
      
        
  }
  fdr.values <- p.adjust(intersection.table[,2], method = "BH")
  intersection.table[,3] <- fdr.values
  write.table(intersection.table, 
              file=paste0("top_parameters_vs_clusters/intersections_", names(top.genes[i]), as.character(threshold),".txt"), 
              sep="\t", row.names = FALSE, quote = FALSE)
}


# degree.intersections <- read.table(file="topvalues_clusters/intersections_Degree0.9.txt", header = TRUE, sep = "\t")
# head(degree.intersections)
# degree.intersections$fdr
# degree.intersections$p.value
# 
# degree.intersections <- read.table(file="topvalues_clusters/intersections_Degree0.7.txt", header = TRUE, sep = "\t")
# head(degree.intersections)
# degree.intersections$fdr
# degree.intersections$p.value





##Loop to intersect the clusters of genes and output of network motifs (FFL)####

#Get the output of the motifs
ffl.data <- read.table(file = "../../data/motifs/feedforward_loops_with_multiple_output.tsv", 
                       header = TRUE, sep = "\t")
head(ffl.data)
ffl.outputs <- ffl.data$output

feedback.data <- read.table(file = "../../data/motifs/feedback_loops_with_multiple_output.tsv", 
                       header = TRUE, sep = "\t")
head(feedback.data)
feedback.outputs <- feedback.data$output



#Initialize matrix to store the results
intersection.table <- matrix(ncol=6, nrow = length(ffl.outputs))
colnames(intersection.table) <- c("master_regulator 1","master_regulator_2", "p-value", "fdr", "enrichment", "Intersection Genes") 
head(intersection.table)
i <- 1
j <- 2    
for (i in 1:length(brc1.clusters))
{
  for (j in 1:length(ffl.outputs))
  {
    
    
    current.cluster <- brc1.clusters[[i]]
    current.output <- strsplit(x = as.character(ffl.outputs[[j]]), split = ",")[[1]]
    
    
    print(paste0("TEST",j))
    result <- intersect2sets(set1 = current.cluster, set2 = current.output, alias = alias, gene.descriptions = description)
    p.value <- result[1][[1]]
    enrichment <- result[2][[1]]
    intersect.genes <- result[3][[1]]$intersection.genes
    
    
    
    intersection.table[j,1] <- as.character(ffl.data[j,1][1])
    intersection.table[j,2] <- as.character(ffl.data[j,2][1])
    intersection.table[j,3] <- p.value
    intersection.table[j,5] <- enrichment
    intersection.table[j,6] <- paste(intersect.genes, collapse= ",")
    
    
  }
  fdr.values <- p.adjust(intersection.table[,3], method = "BH")
  intersection.table[,4] <- fdr.values
  write.table(intersection.table, 
              file=paste0("ffl_outputs_vs_clusters/ffl_outputs_vs_cluster_", i,".txt"), 
              sep="\t", row.names = FALSE, quote = FALSE)
}


######Intersection between BRC1 target genes and all the TFs target genes#########
tf.names <- colnames(brc1.network.data)[21:57]
brc1.targets <- gene.names[brc1.network.data$BRC1 == 1]


#Initialize matrix to store the results
intersection.table <- matrix(ncol=5, nrow = length(tf.names))
colnames(intersection.table) <- c("Transcription factor", "p-value", "fdr", "enrichment", "Intersection Genes") 
head(intersection.table)
   
for (i in 1:length(tf.names))
{
    current.target.genes <- gene.names[brc1.network.data[[i+20]] == 1]
    
    result <- intersect2sets(set1 = brc1.targets, set2 = current.target.genes, alias = alias, gene.descriptions = description)
    p.value <- result[1][[1]]
    enrichment <- result[2][[1]]
    intersect.genes <- result[3][[1]]$intersection.genes
    
    
    
    intersection.table[i,1] <- colnames(brc1.network.data)[i+20]
    intersection.table[i,2] <- p.value
    intersection.table[i,4] <- enrichment
    intersection.table[i,5] <- paste(intersect.genes, collapse= ",")
    
    
}
  
fdr.values <- p.adjust(intersection.table[,2], method = "BH")
intersection.table[,3] <- fdr.values
write.table(intersection.table, 
            file="BRC1_targets_vs_all_TFs_targets.txt", 
            sep="\t", row.names = FALSE, quote = FALSE)


######Venn diagrams to visualize intersection between target genes####

library(VennDiagram)

intersect.tfs <- intersect(brc1.targets, current.target.genes)

grid.newpage()
draw.pairwise.venn(area1 = length(brc1.targets),
                   area2 = length(current.target.genes),
                   cross.area = length(intersect.tfs),
                   category = c("BRC1 targets", "Selected TF targets"),
                   fill=c("mediumseagreen","darkred"), 
                   cat.pos = c(0,190), cat.dist = rep(0.05, 1.3), cat.cex=1.6,
                   lwd=3,cex=1.4, fontfamily="arial",
                   cat.fontfamily="arial",fontface="bold",
                   cat.fontface="bold")

  