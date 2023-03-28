## Load igraph library
library(igraph)

##Load the network data
network.data <- read.table(file="brc1_network.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "",comment.char = "%")
rownames(network.data) <- network.data$names

tfs.names <- c("ATAF1", "bZIP52", "PIF3", "MYB3", "ZAT10", 
               "ERF055", "VIP1", "ERF014", "NAC018", "NAP", 
               "ANAC032", "RGA", "HB21", "HB6", "SVP", 
               "ABI5", "IBH1", "BRC1", "ERF035", "GBF2", 
               "PDF2", "HSFB2B", "TCX2", "WRKY18", "ABF3", 
               "HB40", "ZAT6", "DREB2A", "MYB56", "ERF003", 
               "NAC6.ORE1", "HAT2", "SPCH", "DOF5_4", "NAC102", 
               "HB53", "HEC1","GBF3", "HSFB2A")

tf.ids <- c("AT1G01720", "AT1G06850" ,"AT1G09530", "AT1G22640", "AT1G27730",
            "AT1G36060", "AT1G43700", "AT1G44830", "AT1G52880", "AT1G69490",
            "AT1G77450", "AT2G01570", "AT2G18550", "AT2G22430", "AT2G22540",
            "AT2G36270", "AT2G43060", "AT3G18550", "AT3G60490", "AT4G01120",
            "AT4G04890", "AT4G11660", "AT4G14770", "AT4G31800", "AT4G34000",
            "AT4G36740", "AT5G04340", "AT5G05410", "AT5G17800", "AT5G25190",
            "AT5G39610", "AT5G47370", "AT5G53210", "AT5G60850", "AT5G63790",
            "AT5G66700", "AT5G67060", "AT2G46270", "AT5G62020")

names(tfs.names) <- tf.ids
names(tf.ids) <- tfs.names

tf.ids["BRC1"]

setdiff(tfs.names,colnames(network.data))
setdiff(colnames(network.data),tfs.names)

## Generate adjacency matrix
tf.targets <- network.data[,tfs.names]
head(tf.targets)
colnames(tf.targets) <- tf.ids[colnames(tf.targets)]
tf.targets <- as.matrix(tf.targets)
is.matrix(tf.targets)

network.genes <- rownames(network.data)
length(network.genes)

adjacency.matrix <- matrix(0,nrow=length(network.genes),ncol=length(network.genes))
dim(adjacency.matrix)
rownames(adjacency.matrix) <- network.genes
colnames(adjacency.matrix) <- network.genes

dim(adjacency.matrix[,colnames(tf.targets)])
dim(tf.targets)
adjacency.matrix[,colnames(tf.targets)] <- tf.targets
dim(adjacency.matrix)
adjacency.matrix <- t(adjacency.matrix)

  
brc1.graph <- graph.adjacency(adjmatrix = adjacency.matrix,mode = "directed")
vertex.names <- V(brc1.graph)$name
number.nodes <- length(vertex.names)
out.degree <- degree(brc1.graph,mode = "out")
number.tfs <- sum(out.degree != 0)
tfs.out.degree <- out.degree[out.degree != 0]

in.degree <- degree(brc1.graph,mode = "in")

sum(table(in.degree))

hist(in.degree)
hist(degree(brc1.graph),breaks=500,xlim=c(0,50))

write.graph(graph = brc1.graph,file = "brc1net.gml",format = "gml")

## Network motif with three nodes

## graph.motifs inputs consists of a network and a subgraph size k
## (only k= 3 or 4 are supported). It outputs the number of occurrences
## of any subgraph of size k in the given network.

occurrency.subgraph.three.nodes.in.brc1 <- graph.motifs(brc1.graph, size=3)
occurrency.subgraph.three.nodes.in.brc1
length(occurrency.subgraph.three.nodes.in.brc1)
write(occurrency.subgraph.three.nodes.in.brc1,file="motifs/occurency_subgraph_three_nodes_in_brc1.txt",ncolumns = 1)

## Subgraphs of size 3
plot.igraph(graph.isocreate(size=3, number=0))
plot.igraph(graph.isocreate(size=3, number=1))
plot.igraph(graph.isocreate(size=3, number=2))
plot.igraph(graph.isocreate(size=3, number=3))
plot.igraph(graph.isocreate(size=3, number=4))
plot.igraph(graph.isocreate(size=3, number=5))
plot.igraph(graph.isocreate(size=3, number=6))
plot.igraph(graph.isocreate(size=3, number=7))
plot.igraph(graph.isocreate(size=3, number=8))
plot.igraph(graph.isocreate(size=3, number=9))
plot.igraph(graph.isocreate(size=3, number=10))
plot.igraph(graph.isocreate(size=3, number=11))
plot.igraph(graph.isocreate(size=3, number=12))
plot.igraph(graph.isocreate(size=3, number=13))
plot.igraph(graph.isocreate(size=3, number=14))
plot.igraph(graph.isocreate(size=3, number=15))

## Generate randomisation
number.randomisation <- 10000

random.autorregulations <- vector(mode = "numeric", length = number.randomisation)
motifs.3.random.graph <- matrix(0,nrow=number.randomisation, ncol=16)
moments.random.graph <- matrix(0,nrow=number.randomisation,ncol=5)


random.tfs <- sample(x=1:number.nodes, size=number.tfs,replace=FALSE)

install.packages("moments")
library(moments)

gsize(brc1.graph)
gorder(brc1.graph)
brc1.mean <- mean(degree(brc1.graph))
brc1.var <- var(degree(brc1.graph))
brc1.skewness <- skewness(degree(brc1.graph))

for (i in 1:number.randomisation)
{
 ## Print to keep track of the process 
 print(paste0("motif 3 ",i))
 
 ## Initialize random adjancency matrix with zeros
 current.random.adjacency <- matrix(0,nrow=number.nodes,ncol=number.nodes)
 
 ## Randomly select the nodes that will act as TFs with the same number of 
 ## targets as the original BRC1NET
 random.tfs <- sample(x=1:number.nodes, size=number.tfs,replace=FALSE)
 
 
  for(j in 1:length(random.tfs))
  {
    current.random.adjacency[random.tfs[j], sample(x=1:number.nodes, tfs.out.degree[j], replace=FALSE)] <- 1
  }
  
  current.random.network <- graph.adjacency(adjmatrix = current.random.adjacency,mode = "directed")
  moments.random.graph[i,1] <- mean(degree(current.random.network)) 
  moments.random.graph[i,2] <- var(degree(current.random.network))
  moments.random.graph[i,3] <- 100*abs(var(degree(current.random.network)) - brc1.var)/brc1.var
  moments.random.graph[i,4] <- skewness(degree(current.random.network))
  moments.random.graph[i,5] <- 100*abs(skewness(degree(current.random.network)) - brc1.skewness)/brc1.skewness
  
  #random.autorregulations[i] <- sum(diag(current.random.adjacency))
  motifs.3.random.graph[i,] <- graph.motifs(current.random.network, size=3)
}

write.table(x = motifs.3.random.graph,file = "motifs/motifs_three_random_graph.txt",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(x = moments.random.graph,file = "motifs/random_networks_moments.txt",quote = F,row.names = F,col.names = F,sep = "\t")

random.network.moments <- read.table(file = "motifs/random_networks_moments.txt",header = T)
head(random.network.moments)
summary(random.network.moments[[3]])
summary(random.network.moments[[5]])

motifs.3.random.graph <- read.table(file="motifs/motifs_three_random_graph.txt",header=F,as.is=T)
head(motifs.3.random.graph)

plot.igraph(graph.isocreate(size=3, number=0))


## Test significance for each motif
estimated.p.values <- vector(mode="numeric", length=16)
for(i in 1:16)
{
  estimated.p.values[i] <- sum(motifs.3.random.graph[,i] > occurrency.subgraph.three.nodes.in.brc1[i])/number.randomisation
}

indeces.significant.motifs.3 <- which(estimated.p.values < 0.01) - 1
write(x = indeces.significant.motifs.3,file = "motifs/indeces_significant_motifs_3.txt",ncolumns = 1)

## Load three nodes motifs indeces
motifs.3.ind <- read.table(file="indeces_significant_motifs_3.txt")[[1]]

i <- 1
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[i]))
occurrency.subgraph.three.nodes.in.brc1[motifs.3.ind[i]+1]
mean(motifs.3.random.graph[,motifs.3.ind[i]+1])
sd(motifs.3.random.graph[,motifs.3.ind[i]+1])
sum(motifs.3.random.graph[,motifs.3.ind[i]+1] < motifs.3.random.graph[,motifs.3.ind[i]+1])/1000
i <- i + 1



plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[2]))

## Load three nodes motifs occurrences in attractor
occurrences.3 <- read.table(file="occurency_subgraph_three_nodes_in_brc1.txt")[[1]]

## Motif 1
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[3]))
occurrences.3[motifs.3.ind[3] + 1]

maps.motif.1 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[3]), 
                                        brc1.graph, all.maps=TRUE)[["maps"]]

maps.motif.1

current.instance <- names(maps.motif.1[[1]])
master.regulator <- tfs.names[current.instance[3]]
output.gene.1 <- tfs.names[current.instance[2]]
output.gene.2 <- tfs.names[current.instance[1]]

output.tfs <- paste(sort(c(output.gene.1, output.gene.2)),collapse=":")

motif.1.instances <- matrix(c(master.regulator,output.tfs),nrow=1)
colnames(motif.1.instances) <- c("master_regulator","outputs")

for(i in 2:length(maps.motif.1))
{
  current.instance <- names(maps.motif.1[[i]])
  master.regulator <- tfs.names[current.instance[3]]
  output.gene.1 <- tfs.names[current.instance[2]]
  output.gene.2 <- tfs.names[current.instance[1]]
  output.tfs <- paste(sort(c(output.gene.1, output.gene.2)),collapse=":")
  
  
  ffl.to.add <- which((motif.1.instances[,"master_regulator"] == master.regulator))

  if(length(ffl.to.add) != 0)
  {
    motif.1.instances[ffl.to.add,2] <- paste(unique(c(
      strsplit(motif.1.instances[ffl.to.add,2],split=",")[[1]],output.tfs)),collapse=",")
  } else
  {
    motif.1.instances <- rbind(motif.1.instances,c(master.regulator,
                                                   output.tfs))
  }
}

head(motif.1.instances)
motif.1.instances[motif.1.instances[,1] == "BRC1",]
length(strsplit(motif.1.instances[motif.1.instances[,1] == "BRC1",2],",")[[1]])
nrow(motif.1.instances)

write.table(x = motif.1.instances,file = "motifs/regulated_feedback_loop_20201016.tsv",quote = F,sep = "\t",row.names = F)

## Motif 2
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[2]))
occurrences.3[motifs.3.ind[2] + 1]

maps.motif.2 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[2]), 
                                        brc1.graph, all.maps=TRUE)[["maps"]]

maps.motif.2

current.instance <- names(maps.motif.2[[1]])
master.regulator <- tfs.names[current.instance[3]]
second.regulator <- tfs.names[current.instance[2]]
output.gene <- current.instance[1]

regulators <- paste(c(master.regulator, second.regulator),collapse=",")

motif.2.instances <- matrix(c(regulators,output.gene),nrow=1)
colnames(motif.2.instances) <- c("regulators","output")

for(i in 2:length(maps.motif.2))
{
  current.instance <- names(maps.motif.2[[i]])
  master.regulator <- tfs.names[current.instance[3]]
  second.regulator <- tfs.names[current.instance[2]]
  output.gene <- current.instance[1]
  
  regulators <- paste(c(master.regulator, second.regulator),collapse=",")
  print(regulators)
  
  ffl.to.add <- which((motif.2.instances[,"regulators"] == regulators))
  
  if(length(ffl.to.add) != 0)
  {
    motif.2.instances[ffl.to.add,2] <- paste(motif.2.instances[ffl.to.add,2],output.gene,sep=",")
  } else
  {
    motif.2.instances <- rbind(motif.2.instances,c(regulators,output.gene))
  }
}

motif.2.instances[,1]
nrow(motif.2.instances)
length(strsplit(motif.1.instances[motif.1.instances[,1] == "BRC1",2],",")[[1]])

write.table(x = motif.2.instances,file = "motifs/feedforward_loop_with_multiple_output.tsv",quote = F,sep = "\t",row.names = F)

## Motif 3
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[7]))
occurrences.3[motifs.3.ind[7] + 1]

maps.motif.3 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[7]), 
                                        brc1.graph, all.maps=TRUE)[["maps"]]

maps.motif.3

current.instance <- names(maps.motif.3[[1]])
master.regulator1 <- tfs.names[current.instance[3]]
master.regulator2 <- tfs.names[current.instance[1]]
output.gene <- current.instance[2]

master.regulators <- paste(sort(c(master.regulator1, master.regulator2)),collapse=",")

motif.3.instances <- matrix(c(master.regulators,output.gene),nrow=1)
colnames(motif.3.instances) <- c("master_regulators","output")

i <- 2
for(i in 2:length(maps.motif.3))
{
  current.instance <- names(maps.motif.3[[i]])
  master.regulator1 <- tfs.names[current.instance[3]]
  master.regulator2 <- tfs.names[current.instance[1]]
  output.gene <- current.instance[2]
  
  master.regulators <- paste(sort(c(master.regulator1, master.regulator2)),collapse=",")
  
  ffl.to.add <- which((motif.3.instances[,"master_regulators"] == master.regulators))
  
  if(length(ffl.to.add) != 0)
  {
    motif.3.instances[ffl.to.add,2] <- paste(motif.3.instances[ffl.to.add,2],output.gene,sep=",")
  } else
  {
    motif.3.instances <- rbind(motif.3.instances,c(master.regulators,output.gene))
  }
}

head(motif.3.instances)

motif.3.instances[,1]
nrow(motif.3.instances)
length(strsplit(motif.1.instances[motif.1.instances[,1] == "BRC1",2],",")[[1]])

write.table(x = motif.3.instances,file = "motifs/feedback_loop_with_multiple_output.tsv",quote = F,sep = "\t",row.names = F)


## Motif 3
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[8]))
occurrences.3[motifs.3.ind[8] + 1]


tf.ids <- setdiff(tf.ids,"AT3G18550" )

## 
i <- 1
for(i in 1:10)
{
  print(i)
  tfs.to.remove <- combn(x = tf.ids, m = i, simplify = FALSE)
  print("done generation")
  
  if(length(tfs.to.remove) > 50000)
  {
    tfs.to.remove <- tfs.to.remove[1:50000]
  }
  
  mean.shortest.paths <- vector(mode = "numeric",length = length(tfs.to.remove))
  number.components <- vector(mode = "numeric",length = length(tfs.to.remove))
  for(j in 1:length(tfs.to.remove))
  {
    new.brc1.network <- delete.vertices(graph = brc1.graph,v = tfs.to.remove[[j]])
    mean.shortest.paths[j] <- mean_distance(graph = new.brc1.network)
    number.components[j] <- count_components(graph = new.brc1.network)
  }
  print("done computation")
  write.table(x = mean.shortest.paths,file = paste(c("results_robustness/brc1_mean_shortest_paths_",i,".txt"),collapse=""),quote = F,row.names = F,col.names = F)
  write.table(x = number.components,file = paste(c("results_robustness/brc1_number_components_",i,".txt"),collapse=""),quote = F,row.names = F,col.names = F)
}

dist.1 <- read.table(file="results_robustness/mean_shortest_paths_1.txt",header = F,as.is = T)[[1]]
dist.2 <- read.table(file="results_robustness/mean_shortest_paths_2.txt",header = F,as.is = T)[[1]]
dist.3 <- read.table(file="results_robustness/mean_shortest_paths_3.txt",header = F,as.is = T)[[1]]
dist.4 <- read.table(file="results_robustness/mean_shortest_paths_4.txt",header = F,as.is = T)[[1]]
dist.5 <- read.table(file="results_robustness/mean_shortest_paths_5.txt",header = F,as.is = T)[[1]]
dist.6 <- read.table(file="results_robustness/mean_shortest_paths_6.txt",header = F,as.is = T)[[1]]
dist.7 <- read.table(file="results_robustness/mean_shortest_paths_7.txt",header = F,as.is = T)[[1]]
dist.8 <- read.table(file="results_robustness/mean_shortest_paths_8.txt",header = F,as.is = T)[[1]]
dist.9 <- read.table(file="results_robustness/mean_shortest_paths_9.txt",header = F,as.is = T)[[1]]
dist.10 <- read.table(file="results_robustness/mean_shortest_paths_10.txt",header = F,as.is = T)[[1]]



boxplot(dist.1,dist.2,dist.3,dist.4,dist.5,dist.6,dist.7,dist.8,dist.9,dist.10,outline=F,ylim=c(1.5,2),col="cornflowerblue",
        main="Shortest Paths Lengths",xlab="Number of removed TFs",ylab="Length",cex.lab=1.5,cex.main=2,names=1:9)

plot(c(mean(dist.1),
       mean(dist.2),
       mean(dist.3),
       mean(dist.4),
       mean(dist.5),
       mean(dist.6),
       mean(dist.7),
       mean(dist.8),
       mean(dist.9)),lwd=3,ylim=c(1.5,2),type="o",main="Average ")

dist.1 <- read.table(file="results_robustness/number_components_1.txt",header = F,as.is = T)[[1]]
dist.2 <- read.table(file="results_robustness/number_components_2.txt",header = F,as.is = T)[[1]]
dist.3 <- read.table(file="results_robustness/number_components_3.txt",header = F,as.is = T)[[1]]
dist.4 <- read.table(file="results_robustness/number_components_4.txt",header = F,as.is = T)[[1]]
dist.5 <- read.table(file="results_robustness/number_components_5.txt",header = F,as.is = T)[[1]]
dist.6 <- read.table(file="results_robustness/number_components_6.txt",header = F,as.is = T)[[1]]
dist.7 <- read.table(file="results_robustness/number_components_7.txt",header = F,as.is = T)[[1]]
dist.8 <- read.table(file="results_robustness/number_components_8.txt",header = F,as.is = T)[[1]]
dist.9 <- read.table(file="results_robustness/number_components_9.txt",header = F,as.is = T)[[1]]

boxplot(dist.1,dist.2,dist.3,dist.4,dist.5,dist.6,dist.7,dist.8,dist.9,outline=F)

plot(c(mean(dist.1),
       mean(dist.2),
       mean(dist.3),
       mean(dist.4),
       mean(dist.5),
       mean(dist.6),
       mean(dist.7),
       mean(dist.8),
       mean(dist.9)),lwd=3,ylim=c(0,10),type="o")



dist.1 <- read.table(file="results_robustness/brc1_mean_shortest_paths_1.txt",header = F,as.is = T)[[1]]
dist.2 <- read.table(file="results_robustness/brc1_mean_shortest_paths_2.txt",header = F,as.is = T)[[1]]
dist.3 <- read.table(file="results_robustness/brc1_mean_shortest_paths_3.txt",header = F,as.is = T)[[1]]
dist.4 <- read.table(file="results_robustness/brc1_mean_shortest_paths_4.txt",header = F,as.is = T)[[1]]
dist.5 <- read.table(file="results_robustness/brc1_mean_shortest_paths_5.txt",header = F,as.is = T)[[1]]
dist.6 <- read.table(file="results_robustness/brc1_mean_shortest_paths_6.txt",header = F,as.is = T)[[1]]
dist.7 <- read.table(file="results_robustness/brc1_mean_shortest_paths_7.txt",header = F,as.is = T)[[1]]
dist.8 <- read.table(file="results_robustness/brc1_mean_shortest_paths_8.txt",header = F,as.is = T)[[1]]
dist.9 <- read.table(file="results_robustness/brc1_mean_shortest_paths_9.txt",header = F,as.is = T)[[1]]
dist.10 <- read.table(file="results_robustness/brc1_mean_shortest_paths_10.txt",header = F,as.is = T)[[1]]



boxplot(dist.1,dist.2,dist.3,dist.4,dist.5,dist.6,dist.7,dist.8,dist.9,dist.10,outline=F,ylim=c(1.5,2),col="cornflowerblue",
        main="Shortest Paths Lengths",xlab="Number of removed TFs",ylab="Length",cex.lab=1.5,cex.main=2,names=1:10)


dist.1 <- read.table(file="results_robustness/brc1_number_components_1.txt",header = F,as.is = T)[[1]]
dist.2 <- read.table(file="results_robustness/brc1_number_components_2.txt",header = F,as.is = T)[[1]]
dist.3 <- read.table(file="results_robustness/brc1_number_components_3.txt",header = F,as.is = T)[[1]]
dist.4 <- read.table(file="results_robustness/brc1_number_components_4.txt",header = F,as.is = T)[[1]]
dist.5 <- read.table(file="results_robustness/brc1_number_components_5.txt",header = F,as.is = T)[[1]]
dist.6 <- read.table(file="results_robustness/brc1_number_components_6.txt",header = F,as.is = T)[[1]]
dist.7 <- read.table(file="results_robustness/brc1_number_components_7.txt",header = F,as.is = T)[[1]]
dist.8 <- read.table(file="results_robustness/brc1_number_components_8.txt",header = F,as.is = T)[[1]]
dist.9 <- read.table(file="results_robustness/brc1_number_components_9.txt",header = F,as.is = T)[[1]]
dist.10 <- read.table(file="results_robustness/brc1_number_components_10.txt",header = F,as.is = T)[[1]]


boxplot(dist.1,dist.2,dist.3,dist.4,dist.5,dist.6,dist.7,dist.8,dist.9,dist.10,outline=F,col="salmon",
        main="Conected Components",xlab="Number of removed TFs",ylab="Number of Components",cex.lab=1.5,cex.main=2,names=1:10)



boxplot(dist.1,dist.2,dist.3,dist.4,dist.5,outline=F)

plot(c(mean(dist.1),
       mean(dist.2),
       mean(dist.3),
       mean(dist.4),
       mean(dist.5),
       mean(dist.6),
       mean(dist.7),
       mean(dist.8),
       mean(dist.9),
       mean(dist.10)),lwd=3,ylim=c(0,12),type="o")




library(rje)
powerSet(1:3,m = 2)
help(powerSet)
combn(tf.ids, 2, simplify = FALSE)
