library(igraph)

##Load the network data
network.data <- read.table(file="brc1_network.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "",comment.char = "%")
rownames(network.data) <- network.data$names

tfs.names <- c("NAC002", "AtbZIP52", "PIF3", "MYB3", "ZAT10", 
               "ERF055", "VIP1", "ERF014", "NAC018", "NAP", 
               "ANAC032", "RGA", "ATHB_21", "ATHB_6", "SVP", 
               "ABI5", "IBH1", "BRC1", "ERF035", "GBF2", 
               "PDF2", "HSFB2B", "TCX2", "WRKY18", "ABF3", 
               "ATHB_40", "ZAT6", "DREB2A", "AtMYB56", "ERF003", 
               "NAC6_ORE1", "HAT2", "SPCH", "DOF5_4", "NAC102", 
               "ATHB_53", "HEC1","GBF3")

tf.ids <- c("AT1G01720", "AT1G06850" ,"AT1G09530", "AT1G22640", "AT1G27730",
            "AT1G36060", "AT1G43700", "AT1G44830", "AT1G52880", "AT1G69490",
            "AT1G77450", "AT2G01570", "AT2G18550", "AT2G22430", "AT2G22540",
            "AT2G36270", "AT2G43060", "AT3G18550", "AT3G60490", "AT4G01120",
            "AT4G04890", "AT4G11660", "AT4G14770", "AT4G31800", "AT4G34000",
            "AT4G36740", "AT5G04340", "AT5G05410", "AT5G17800", "AT5G25190",
            "AT5G39610", "AT5G47370", "AT5G53210", "AT5G60850", "AT5G63790",
            "AT5G66700", "AT5G67060", "AT2G46270")

names(tfs.names) <- tf.ids
names(tf.ids) <- tfs.names

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
motifs.3.random.graph <- read.table(file="motifs_three_random_graph.txt",header=F,as.is=T)
head(motifs.3.random.graph)
write.table(x = random.autorregulations,file = "autorregulation_random_graph.txt",quote = F,row.names = F,col.names = F,sep = "\t")


## Test significance for each motif
estimated.p.values <- vector(mode="numeric", length=16)
for(i in 1:16)
{
  estimated.p.values[i] <- sum(motifs.3.random.graph[,i] > occurrency.subgraph.three.nodes.in.brc1[i])/number.randomisation
}

indeces.significant.motifs.3 <- which(estimated.p.values < 0.01) - 1
write(x = indeces.significant.motifs.3,file = "indeces_significant_motifs_3.txt",ncolumns = 1)

## Load three nodes motifs indeces
motifs.3.ind <- read.table(file="indeces_significant_motifs_3.txt")[[1]]
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[2]))

## Load three nodes motifs occurences in attractor
occurrences.3 <- read.table(file="occurency_subgraph_three_nodes_in_brc1.txt")[[1]]

## Motif 1
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[2]))
occurrences.3[motifs.3.ind[2] + 1]

maps.motif.1 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[2]), 
                                        brc1.graph, all.maps=TRUE)[["maps"]]

maps.motif.1

current.instance <- names(maps.motif.1[[1]])
master.regulator <- tfs.names[current.instance[3]]
output.gene.1 <- tfs.names[current.instance[2]]
output.gene.2 <- tfs.names[current.instance[1]]

motif.1.instances <- matrix(c(master.regulator,paste(output.gene.1,output.gene.2,sep=",")),nrow=1)
colnames(motif.1.instances) <- c("master_regulator","outputs")

for(i in 2:length(maps.motif.1))
{
  current.instance <- names(maps.motif.1[[i]])
  master.regulator <- tfs.names[current.instance[3]]
  output.gene.1 <- tfs.names[current.instance[2]]
  output.gene.2 <- tfs.names[current.instance[1]]
  
  ffl.to.add <- which((motif.1.instances[,"master_regulator"] == master.regulator))

  if(length(ffl.to.add) != 0)
  {
    motif.1.instances[ffl.to.add,2] <- paste(unique(c(
      strsplit(motif.1.instances[ffl.to.add,2],split=",")[[1]],output.gene.1,output.gene.2)),collapse=",")
  } else
  {
    motif.1.instances <- rbind(motif.1.instances,c(master.regulator,
                                                   paste(output.gene.1,output.gene.2,sep=",")))
  }
}

head(motif.1.instances)
motif.1.instances[motif.1.instances[,1] == "BRC1",]
length(strsplit(motif.1.instances[motif.1.instances[,1] == "BRC1",2],",")[[1]])

write.table(x = motif.1.instances,file = "regulated_feedback_loop.tsv",quote = F,sep = "\t",row.names = F)

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

write.table(x = motif.2.instances,file = "feedforward_loop_with_multiple_output.tsv",quote = F,sep = "\t",row.names = F)





## Motif 3
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[6]))
occurrences.3[motifs.3.ind[6] + 1]

maps.motif.3 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[6]), 
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

write.table(x = motif.3.instances,file = "feedback_loop_with_multiple_output.tsv",quote = F,sep = "\t",row.names = F)


## Motif 3
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[8]))
occurrences.3[motifs.3.ind[8] + 1]




# No muy interesante

