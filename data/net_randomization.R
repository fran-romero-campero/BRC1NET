############################
#  net_randomization.R
############################
# This script ...
#
#
# author = apascualgarcia.github.io
# 

## Load igraph library
library(igraph)
library(ggplot2)

# Set environment and load local functions ------
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
source('curveball.versions.R')
source('perturbation.R')


# Read file --------

brc1.graph = read.graph(file = "brc1net.gml",format = "gml")
bcr1.adjlist <- get.adjlist(brc1.graph, mode = "out")
brc1.adj = get.adjacency(brc1.graph)

is.simple(brc1.graph) # it includes loops


# Randomize directed network 
Nrnd = c(50, 100, 500, 1000, 5000, 10000, 100000, 200000)
nn = Nrnd[8]
sample <- curveball.randomise.dir(bcr1.adjlist, nn, T, T)
brc1.rnd.adjlist <- graph.adjlist(sample[[nn]])
brc1.rnd.adj <- get.adjacency(brc1.rnd.adjlist)

plot(degree(brc1.graph, mode = "in"), degree(brc1.graph.rnd, mode = "in"))
plot(degree(brc1.graph, mode = "out"), degree(brc1.graph.rnd, mode = "out"))
prod(degree(brc1.graph, mode = "in") == degree(brc1.graph.rnd, mode = "in"))
prod(degree(brc1.graph, mode = "out") == degree(brc1.graph.rnd, mode = "out"))

# compute the proportion of links that are different between observed and random
interval = 100
prop.diff = vector("numeric",length = (nn/interval))
Nedges = sum(brc1.adj) 
index = seq(0,nn,interval); index[1] = 1
for(i in index){ 
  brc1.rnd.adjlist.tmp <- graph.adjlist(sample[[i]])
  brc1.rnd.adj.tmp <- get.adjacency(brc1.rnd.adjlist.tmp)
  prop.diff[i] = sum(brc1.adj != brc1.rnd.adj.tmp)/Nedges
}

# Plot results
df = data.frame(step = index, distance = prop.diff[index])

file.Out = "motifs/Plot_curveball_convergence.pdf"
pdf(file.Out)
gg = ggplot(df, aes(step, distance))+
  geom_point()+theme_bw()
print(gg)
dev.off()


plot(prop.diff)

#ps = perturbation.score(sample,nn)
