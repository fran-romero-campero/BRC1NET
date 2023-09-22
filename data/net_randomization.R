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


# Set environment and load local functions ------
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
source('curveball.versions.R')
source('perturbation.R')


# Read file --------

brc1.graph = read.graph(file = "brc1net.gml",format = "gml")
I <- get.adjlist(brc1.graph, mode = "out")
is.simple(brc1.graph) # it includes loops


# Randomize directed network 
Nrnd = c(50, 100, 500, 1000, 5000, 10000)
nn = Nrnd[6]
sample <- curveball.randomise.dir(I, nn, T, T)
brc1.graph.rnd <- graph.adjlist(sample[[nn]])
plot(degree(brc1.graph, mode = "in"), degree(brc1.graph.rnd, mode = "in"))
plot(degree(brc1.graph, mode = "out"), degree(brc1.graph.rnd, mode = "out"))
prod(degree(brc1.graph, mode = "in") == degree(brc1.graph.rnd, mode = "in"))
prod(degree(brc1.graph, mode = "out") == degree(brc1.graph.rnd, mode = "out"))


ps = perturbation.score(sample,nn)
