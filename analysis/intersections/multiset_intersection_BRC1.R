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
  results <- supertest(x = sets, n = 5778)
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
  results <- supertest(x = sets, n = 5778)
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
intersection.table <- matrix(ncol=6, nrow = length(feedback.outputs))
colnames(intersection.table) <- c("master_regulator 1","master_regulator_2", "p-value", "fdr", "enrichment", "Intersection Genes") 
head(intersection.table)
i <- 1
j <- 2    
for (i in 1:length(brc1.clusters))
{
  for (j in 1:length(feedback.outputs))
  {
    
    
    current.cluster <- brc1.clusters[[i]]
    current.output <- strsplit(x = as.character(feedback.outputs[[j]]), split = ",")[[1]]
    
    
    print(paste0("TEST",j))
    result <- intersect2sets(set1 = current.cluster, set2 = current.output, alias = alias, gene.descriptions = description)
    p.value <- result[1][[1]]
    enrichment <- result[2][[1]]
    intersect.genes <- result[3][[1]]$intersection.genes
    
    
    
    intersection.table[j,1] <- as.character(feedback.data[j,1][1])
    intersection.table[j,2] <- as.character(feedback.data[j,2][1])
    intersection.table[j,3] <- p.value
    intersection.table[j,5] <- enrichment
    intersection.table[j,6] <- paste(intersect.genes, collapse= ",")
    
    
  }
  fdr.values <- p.adjust(intersection.table[,3], method = "BH")
  intersection.table[,4] <- fdr.values
  write.table(intersection.table, 
              file=paste0("feedback_outputs_vs_clusters/feedback_outputs_vs_cluster_", i,".txt"), 
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




#####Intersections between binding regions in DNA (BED files)####
#Reading the bed files of the transcription factors
peaks1 <- read.table(file = "../../../web_apps/peak_visualizer/data/bed_files/ELF3_ZT0_1_peaks.narrowPeak")
head(peaks1)

peaks2 <- read.table(file = "../../../web_apps/peak_visualizer/data/bed_files/PIF3_peaks.narrowPeak")
head(peaks2)

peaks3 <- read.table(file = "../../../web_apps/peak_visualizer/data/bed_files/PRR9_1_peaks.narrowPeak")
head(peaks3)

peaks.list <- list(peaks1, peaks2, peaks3)

length.sets <- sapply(X = peaks.list, FUN = nrow)


peaks.set1 <- peaks1
peaks.set2 <- peaks2
peaks.set2 <- random.peaks2

#intersectBed functions
intersectBed.1 <- function(peaks.set1, peaks.set2)
{
  intersection <- matrix(ncol = 3, nrow=0 )
  current.intersection <- matrix(ncol = 3 )
  for (i in 1:nrow(peaks.set1))
  {
    #Set the current peak values of set1
    current.chr <- as.numeric(peaks.set1[i,1])
    current.start <- peaks.set1[i,2]
    current.end <- peaks.set1[i,3]
    #Checking if there is intersection between the current peak and any peak of set2
    option1 <- nrow(subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.start & peaks.set2[,3]>=current.start))
    option2 <- nrow(subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.end & peaks.set2[,3]>=current.end))
    
    # print(i)
    
    if(option1+option2 > 0)
    {
      # print("HIT")
      if(option1>0)
      {
        hit.peak2 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.start & peaks.set2[,3]>=current.start)
        current.intersection[1,1] <- current.chr
        current.intersection[1,2] <- current.start
        current.intersection[1,3] <- hit.peak2[1,3]
        
      }else
      {
        hit.peak2 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.end & peaks.set2[,3]>=current.end)
        current.intersection[1,1] <- current.chr
        current.intersection[1,2] <- hit.peak2[1,2]
        current.intersection[1,3] <- current.end
      }
      
      intersection <- rbind(intersection, current.intersection)
    }
  }
  return(intersection)
}

intersectBed <- function(peaks.set1, peaks.set2)
{
  intersection <- matrix(ncol = 3, nrow=min(nrow(peaks.set1),nrow(peaks.set2)) )
  j <- 0
  for (i in 1:nrow(peaks.set1))
  {
    #Set the current peak values of set1
    current.chr <- as.numeric(peaks.set1[i,1])
    current.start <- peaks.set1[i,2]
    current.end <- peaks.set1[i,3]
    #Checking if there is intersection between the current peak and any peak of set2
    hit.1 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.start & peaks.set2[,3]>=current.start)
    option1 <- nrow(hit.1)
    
    if(option1>0)
    {
      j <- j + 1
      intersection[j,1] <- current.chr
      intersection[j,2] <- current.start
      intersection[j,3] <- hit.1[1,3]
    } else 
    {
      hit.2 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.end & peaks.set2[,3]>=current.end)    
      option2 <- nrow(hit.2)
      
      if(option2 > 0)
      {
        j <- j + 1
        intersection[j,1] <- current.chr
        intersection[j,2] <- hit.2[1,2]
        intersection[j,3] <- current.end
      }
    }
  }
  return(list(intersection=intersection[1:j,],size=j))
}

#Function to randomize peaks
randomize.peaks <- function(input.peaks,chr.lengths)
{
  random.numbers <- runif(n = nrow(input.peaks))
  random.starts <- ceiling(random.numbers * chr.lengths[input.peaks$V1])
  random.ends <- random.starts + (input.peaks$V3 - input.peaks$V2)
  
  random.peaks <- matrix(c(input.peaks$V1,random.starts,random.ends),ncol=3) 
  return(random.peaks)
}



# The intersectBed function allow to get the intersections between two bed files. If you want to perform the
# the intersection between three bed files, you can do it in a consecutive manner. 

first <- intersectBed(peaks1, peaks2)[[1]]
nrow(first)
second <- intersectBed(first, peaks3)
nrow(second)


## Permutation of peaks.set2 (random.peaks2) and comparing with peaks.set1 ####
chromosomes.length <- read.table(file="../../../web_apps/peak_visualizer/data/bed_files/atha_chr_lengths.txt",as.is=T)[[1]]
number.randomisation <- 20 #100000

random.intersections <- vector(mode = "numeric",length=number.randomisation) #Creating vector
for(j in 1:number.randomisation)
{
  # print(j)
  # random.peaks2 <- matrix(nrow=nrow(peaks2),ncol=3) #Matriz con 3 columnas, una para el cromosoma, otra para el comienzo y otra para el final de la marca aleatoria.
  # for(i in 1:nrow(peaks2))
  # {
  #   current.chr <- peaks2[i,1][[1]] #Chr de la iésima marca real
  #   current.start <- peaks2[i,2] #Start de la iésima marca real
  #   current.end <- peaks2[i,3] #End de la iésima marca real
  #   current.length <- current.end - current.start #Longitud de la iésima marca real
  #   
  #   chr.length <- chromosomes.length[current.chr] #Length del actual cromosoma
  #   #Ahora genero los mismos datos para marcas aleatorias
  #   random.start <- floor(runif(n = 1,min = 1,max = chr.length))
  #   random.end <- random.start + current.length
  #   
  #   random.peaks2[i,1] <- current.chr
  #   random.peaks2[i,2] <- random.start
  #   random.peaks2[i,3] <- random.end
  # }
  
  random.peaks2 <- randomize.peaks(input.peaks = peaks2, chr.lengths = chromosomes.length)
  random.intersections[j] <- intersectBed(peaks.set1 = peaks1, peaks.set2 = random.peaks2)[[2]]
  
  p.value <- sum(random.intersections > nrow(first)) / number.randomisation
  p.value 

}


##Loop to check the intersection of binding regions (bed files) between all the transcription factors together and store the results in a table####
chromosomes.length <- read.table(file="../../../web_apps/peak_visualizer/data/bed_files/atha_chr_lengths.txt",as.is=T)[[1]]
number.randomisation <- 1000
bed.files <- list.files(path = "../../../web_apps/peak_visualizer/data/bed_files/", pattern = "peaks.narrowPeak")

combinations <- expand.grid(bed.files, bed.files)
bed.intersections <- matrix(ncol = 6, nrow = nrow(combinations))
colnames(bed.intersections) <- c("TF1", "TF2", "p-value", "fdr", "number of intersections", "Genes" )


txdb <- TxDb.Athaliana.BioMart.plantsmart28
i <- 76
total.tests <- nrow(combinations)
# total.tests <- 50

# Start the clock!
ptm <- proc.time()

for (i in 1:total.tests)
# for (i in 1:nrow(combinations))
{
  # print(paste0("test number ", i, " of ", nrow(combinations)))
  print(paste0((i/total.tests)*100, " %"))
  peaks1 <- read.table(file = paste0("../../../web_apps/peak_visualizer/data/bed_files/", combinations[i,1]))
  peaks2 <- read.table(file = paste0("../../../web_apps/peak_visualizer/data/bed_files/", combinations[i,2]))
  real.intersection <- intersectBed(peaks.set1 = peaks1, peaks.set2 = peaks2)
  if (real.intersection[[2]] > 0)
  {
    random.intersections <- vector(mode = "numeric",length=number.randomisation) #Creating vector
    for(j in 1:number.randomisation)
    {
      print(j)
      # random.peaks2 <- matrix(nrow=nrow(peaks2),ncol=3) #Matriz con 3 columnas, una para el cromosoma, otra para el comienzo y otra para el final de la región aleatoria.
      # for(k in 1:nrow(peaks2))
      # {
      #   current.chr <- peaks2[k,1][[1]] #Chr de la iésima marca real
      #   current.start <- peaks2[k,2] #Start de la iésima marca real
      #   current.end <- peaks2[k,3] #End de la iésima marca real
      #   current.length <- current.end - current.start #Longitud de la iésima marca real
      #   
      #   chr.length <- chromosomes.length[current.chr] #Length del actual cromosoma
      #   #Ahora genero los mismos datos para regiones aleatorias
      #   random.start <- floor(runif(n = 1,min = 1,max = chr.length))
      #   random.end <- random.start + current.length
      #   
      #   random.peaks2[k,1] <- current.chr
      #   random.peaks2[k,2] <- random.start
      #   random.peaks2[k,3] <- random.end
      # }
      
      random.peaks2 <- randomize.peaks(input.peaks = peaks2, chr.lengths = chromosomes.length)
      random.intersections[j] <- intersectBed(peaks.set1 = peaks1, peaks.set2 = random.peaks2 )[[2]] 
      
      
    }
    
    p.value <- sum(random.intersections > real.intersection[[2]]) / number.randomisation
    if( p.value == 0)
    {
      p.value <- 1/number.randomisation
    }
    
    if (real.intersection[[2]]==1)
    {
      real.intersection[[1]] <- matrix(real.intersection[[1]], ncol = 3)
    }
    
    colnames(real.intersection[[1]]) <- c("chromosome", "start", "end")
    granges.intersection <- makeGRangesFromDataFrame(real.intersection[[1]],
                                                     keep.extra.columns=FALSE,
                                                     ignore.strand=FALSE,
                                                     seqinfo=NULL,
                                                     seqnames.field="chromosome",
                                                     start.field="start",
                                                     end.field="end",
                                                     starts.in.df.are.0based=FALSE)
    
    
    
    peakAnno <- annotatePeak(granges.intersection, tssRegion=c(-2000, 0),
                             TxDb=txdb, annoDb="org.At.tair.db")
    
    annot.peaks <- as.data.frame(peakAnno)
    # target.genes <- subset(annot.peaks, distanceToTSS >= 2000 | distanceToTSS <= -2000)$geneId
    target.genes <- subset(annot.peaks, annotation != "Downstream (<1kb)" & annotation != "Downstream (1-2kb)"
                           & annotation != "Distal Intergenic" & annotation != "Downstream (2-3kb)")$geneId
    
    target.genes <- paste(target.genes, collapse = ",")
    
    bed.intersections[i,1] <- strsplit(x = as.character(combinations[i,1]), split = "_peaks")[[1]][1]
    bed.intersections[i,2] <- strsplit(x = as.character(combinations[i,2]), split = "_peaks")[[1]][1]
    bed.intersections[i,3] <- p.value
    bed.intersections[i,5] <- nrow(real.intersection)
    bed.intersections[i,6] <- target.genes
    
  } else 
  {
    bed.intersections[i,1] <- strsplit(x = as.character(combinations[i,1]), split = "_peaks")[[1]][1]
    bed.intersections[i,2] <- strsplit(x = as.character(combinations[i,2]), split = "_peaks")[[1]][1]
    bed.intersections[i,3] <- NA
    bed.intersections[i,5] <- "No intersection"
    bed.intersections[i,6] <- NA
  }
  
  

}

write.table(bed.intersections, file = "bed_intersections.txt", sep = "\t", row.names = FALSE)
# Stop the clock
proc.time() - ptm

library(TxDb.Athaliana.BioMart.plantsmart28)
library(org.At.tair.db)
library(ChIPseeker)

txdb <- TxDb.Athaliana.BioMart.plantsmart28
colnames(real.intersection) <- c("chromosome", "start", "end")
granges.intersection <- makeGRangesFromDataFrame(real.intersection,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field="chromosome",
                         start.field="start",
                         end.field="end",
                         starts.in.df.are.0based=FALSE)



peakAnno <- annotatePeak(granges.intersection, tssRegion=c(-2000, 0),
                         TxDb=txdb, annoDb="org.At.tair.db")

annot.peaks <- as.data.frame(peakAnno)
# head(annot.peaks)
# unique(annot.peaks$annotation)
# target.genes <- subset(annot.peaks, distanceToTSS >= 2000 | distanceToTSS <= -2000)$geneId
target.genes <- subset(annot.peaks, annotation != "Downstream (<1kb)" & annotation != "Downstream (1-2kb)"
                       & annotation != "Distal Intergenic" & annotation != "Downstream (2-3kb)")$geneId
target.genes <- paste(target.genes, collapse = ",")

# cpsets(x = nrow(second), L = length.sets, n = sum(length.sets), 5778, lower.tail = FALSE) ##DUda, cuánto es n (population size)


