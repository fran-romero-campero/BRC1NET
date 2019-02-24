## This script performs an identification of network motifs in BRC1 downstream network

# Authors: Pedro de los Reyes
##         Francisco J. Romero-Campero
# 
# Contact: Francisco J. Romero-Campero - fran@us.es

## Load library and graph
library(igraph)

## Load BRC1 downstream network and extract gene names
brc1.graph <- read.graph(file="../BRC1_transcriptional_network.graphml", format = "graphml")
vertex.names <- V(brc1.graph)$name
number.nodes <- length(vertex.names)

## This network is not scale-free
hist(degree(brc1.graph),breaks = seq(from=0,to=600,by=2),xlim=c(0,100))

hist(degree(brc1.graph,mode = "out"))
hist(degree(brc1.graph,mode = "in"))

## Extract number of tfs and their out degree
out.degree <- degree(brc1.graph,mode = "out")
number.tfs <- sum(out.degree != 0)
write(names(which(out.degree != 0)),file="TFs_agi.txt")

tfs.out.degree <- out.degree[out.degree != 0]

## Network motif with three nodes

## graph.motifs inputs consists of a network and a subgraph size k
## (only k= 3 or 4 are supported). It outputs the number of occurencies
## of any subgraph of size k in the given network.

occurrency.subgraph.three.nodes.in.brc1 <- graph.motifs(brc1.graph, size=3)
occurrency.subgraph.three.nodes.in.brc1
length(occurrency.subgraph.three.nodes.in.brc1)
write(occurrency.subgraph.three.nodes.in.brc1,file="occurency_subgraph_three_nodes_in_brc1.txt",ncolumns = 1)

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
number.randomisation <- 1000

random.autorregulations <- vector(mode = "numeric", length = number.randomisation)
motifs.3.random.graph <- matrix(0,nrow=number.randomisation, ncol=16)

for (i in 1:number.randomisation)
{
  print(paste0("motif 3 ",i))
  random.tfs <- sample(x=1:number.nodes, size=number.tfs,replace=FALSE)
  
  current.random.adjacency <- matrix(0,nrow=number.nodes,ncol=number.nodes)

  for(j in 1:length(random.tfs))
  {
    current.random.adjacency[random.tfs[j], sample(x=1:number.nodes, tfs.out.degree[j], replace=FALSE)] <- 1
  }
  
  current.random.network <- graph.adjacency(adjmatrix = current.random.adjacency,mode = "directed")
  
  random.autorregulations[i] <- sum(diag(current.random.adjacency))
  motifs.3.random.graph[i,] <- graph.motifs(current.random.network, size=3)
}

write.table(x = motifs.3.random.graph,file = "motifs_three_random_graph.txt",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(x = random.autorregulations,file = "autorregulation_random_graph.txt",quote = F,row.names = F,col.names = F,sep = "\t")


## Test significance for each motif
estimated.p.values <- vector(mode="numeric", length=16)
for(i in 1:16)
{
  estimated.p.values[i] <- sum(motifs.3.random.graph[,i] > occurrency.subgraph.three.nodes.in.brc1[i])/number.randomisation
}

indeces.significant.motifs.3 <- which(estimated.p.values < 0.001) - 1
write(x = indeces.significant.motifs.3,file = "indeces_significant_motifs_3.txt",ncolumns = 1)

## Load three nodes motifs indeces
motifs.3.ind <- read.table(file="indeces_significant_motifs_3.txt")[[1]]

## Load three nodes motifs occurences in attractor
occurrences.3 <- read.table(file="occurency_subgraph_three_nodes_in_brc1.txt")[[1]]

## Motif number 1
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[1]))
occurrences.3[motifs.3.ind[1] + 1]
# No muy interesante

## Motif number 2
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[2]))
occurrences.3[motifs.3.ind[2] + 1]
# Feed forward loop

maps.motif.2 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[2]), 
                                        brc1.graph, all.maps=TRUE)[["maps"]]

maps.motif.2

current.instance <- names(maps.motif.2[[1]])
master.regulator <- current.instance[3]
second.regulator <- current.instance[2]
output.gene <- current.instance[1]

motif.2.instances <- matrix(c(master.regulator,second.regulator,output.gene),nrow=1)
colnames(motif.2.instances) <- c("master_regulator","second_regulator","output")

i <- 2
for(i in 2:length(maps.motif.2))
{
  current.instance <- names(maps.motif.2[[i]])
  master.regulator <- current.instance[3]
  second.regulator <- current.instance[2]
  output.gene <- current.instance[1]
  
  ffl.to.add <- which((motif.2.instances[,"master_regulator"] == master.regulator) & (motif.2.instances[,"second_regulator"] == second.regulator))
  length(which(rep(FALSE,3)))
  
  if(length(ffl.to.add) != 0)
  {
    motif.2.instances[ffl.to.add,3] <- paste(motif.2.instances[ffl.to.add,3],output.gene,sep=",")
  } else
  {
    motif.2.instances <- rbind(motif.2.instances,c(master.regulator,second.regulator,output.gene))
  }
}

head(motif.2.instances)

write.table(x = motif.2.instances,file = "feedforward_loops_with_multiple_output.tsv",quote = F,sep = "\t",row.names = F)

## Motif number 3
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[3]))
occurrences.3[motifs.3.ind[3] + 1]

maps.motif.3 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[3]), 
                                        brc1.graph, all.maps=TRUE)[["maps"]]

current.instance <- names(maps.motif.3[[1]])
master.regulator <- current.instance[3]
target.1 <- current.instance[1]
target.2 <- current.instance[2]
targets <- sort(c(target.1,target.2))

motif.3.instances <- matrix(nrow=length(maps.motif.3),ncol=3)
colnames(motif.3.instances) <- c("master_regulator","target1","target2")

for(i in 1:length(maps.motif.3))
{
  current.instance <- names(maps.motif.3[[i]])
  master.regulator <- current.instance[3]
  target.1 <- current.instance[1]
  target.2 <- current.instance[2]
  targets <- sort(c(target.1,target.2))
  
  motif.3.instances[i,] <- c(master.regulator, targets[1], targets[2])
}

motif.3.instances <- motif.3.instances[!duplicated(motif.3.instances),]
head(motif.3.instances)

write.table(x = motif.3.instances,file = "regulated_feedback_loop_output.tsv",quote = F,sep = "\t",row.names = F)

## Regulated feedback loop

## Motif number 4
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[4]))
occurrences.3[motifs.3.ind[4] + 1]
## No muy interesante

## Motif number 5
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[5]))
occurrences.3[motifs.3.ind[5] + 1]
## Conected Feedback loops 

## Motif number 6
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[6]))
occurrences.3[motifs.3.ind[6] + 1]
##????

## Motif number 7
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[7]))
occurrences.3[motifs.3.ind[7] + 1]
## Feedback loop with output

maps.motif.7 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[7]), 
                                        brc1.graph, all.maps=TRUE)[["maps"]]

current.instance <- names(maps.motif.7[[1]])
regulator.1 <- current.instance[3]
regulator.2 <- current.instance[1]
regulators <- sort(c(regulator.1, regulator.2))
output.gene <- current.instance[2]

motif.7.instances <- matrix(c(regulators[1],regulators[2],output.gene),nrow=1)
colnames(motif.7.instances) <- c("regulator_1","regulator_2","output")

for(i in 2:length(maps.motif.7))
{
  current.instance <- names(maps.motif.7[[i]])
  regulators <- sort(c(current.instance[3],current.instance[1]))
  output.gene <- current.instance[2]
  
  fbl.to.add <- which((motif.7.instances[,"regulator_1"] == regulators[1]) & (motif.7.instances[,"regulator_2"] == regulators[2]))

  if(length(fbl.to.add) != 0)
  {
    motif.7.instances[fbl.to.add,3] <- paste(motif.7.instances[fbl.to.add,3],output.gene,sep=",")
  } else
  {
    motif.7.instances <- rbind(motif.7.instances,c(regulators[1],regulators[2],output.gene))
  }
}

head(motif.7.instances)

write.table(x = motif.7.instances,file = "feedback_loops_with_multiple_output.tsv",quote = F,sep = "\t",row.names = F)







## Motif number 8
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[8]))
occurrences.3[motifs.3.ind[8] + 1]
## ??

## Motif number 9
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[9]))
occurrences.3[motifs.3.ind[9] + 1]


