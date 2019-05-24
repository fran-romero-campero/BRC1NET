## Load library and graph
library(igraph)

## Load BRC1 downstream network and extract gene names
brc1.graph <- read.graph(file="../BRC1_transcriptional_network.graphml", format = "graphml")
vertex.names <- V(brc1.graph)$name
number.nodes <- length(vertex.names)

## This network is not scale-free
hist(degree(brc1.graph),breaks = seq(from=0,to=600,by=2),xlim=c(0,100))
power.law.fit(degree.distribution(brc1.graph))


hist(degree(brc1.graph,mode = "out"))
hist(degree(brc1.graph,mode = "in"))

## Extract number of tfs and their out degree
out.degree <- degree(brc1.graph,mode = "out")
number.tfs <- sum(out.degree != 0)
write(names(which(out.degree != 0)),file="TFs_agi.txt")

tfs.out.degree <- out.degree[out.degree != 0]


## Convert to edge list
library(linkcomm)
brc1.network <- as_edgelist(graph = brc1.graph,names = T)
head(brc1.network)
nrow(brc1.network)

## Extract link communities and visualize them
help("getLinkCommunities")
lc <- getLinkCommunities(network = brc1.network , hcmethod = "ward",directed = T)
print(lc)
save(lc,file = "link_communities_directed.rda")
plot(lc, type = "graph", layout = layout.fruchterman.reingold)
plot(lc, type = "graph", layout = "spencer.circle")
plot(lc, type = "graph", layout = "spencer.circle", shownodesin = 3)
plot(lc, type = "members")
plot(lc, type = "summary")
plot(lc, type = "dend")
