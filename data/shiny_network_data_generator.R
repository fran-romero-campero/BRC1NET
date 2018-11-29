# R script for pre-processing 
# Copyright (C) 2018  Francisco J. Romero-Campero, Pedro de los Reyes
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
# Authors: Francisco J. Romero-Campero
#          Pedro de los Reyes Rodríguez
#          
# 
# Contact: Francisco J. Romero-Campero - fran@us.es 
#         
# Date: November 2018

## A good introduction to XML format can be found here:
## https://www.w3schools.com/xml/xml_whatis.asp

## This script uses the R package XML a good primer can be found here:
## https://www.tutorialspoint.com/r/r_xml_files.htm
## http://www.omegahat.net/RSXML/shortIntro.pdf
## http://www.informit.com/articles/article.aspx?p=2215520

## Input parameters
## Transcription factor, file name containing targets and ZT value
input.xgmml.file <- "BRC1_transcriptional_network.xgmml"

## Load the package required to read XML files.
library("XML")
library("methods")

## Parse the xgmml file making sure the attribute name spaces are kept
result <- xmlTreeParse(file = input.xgmml.file, addAttributeNamespaces = TRUE)

## Extract root node
rootNode <- xmlRoot(result)

## Extract nodes
node.elements <- xmlElementsByTagName(el = rootNode,name = "node")

## Store nodes info
number.nodes <- length(node.elements)

node.attributes <- xmlElementsByTagName(el = node.elements[[1]], name = "att")
number.attributes <- length(node.attributes)
attributes.names <- vector(mode = "character",length=length(number.attributes))
for(i in 1:number.attributes)
{
  attributes.names[i] <- xmlAttrs(node.attributes[[i]])[["name"]]
}

attributes.df <- data.frame(matrix(nrow=number.nodes,ncol=length(attributes.names)))
colnames(attributes.df) <- attributes.names

nodes.names <- vector(mode = "character",length = number.nodes)
nodes.x.pos <- vector(mode = "character",length = number.nodes)
nodes.y.pos <- vector(mode = "character",length = number.nodes)
nodes.color <- vector(mode = "character",length = number.nodes)

for(i in 1:number.nodes)
{
  current.node <- node.elements[[i]]
  nodes.names[i] <- xmlAttrs(current.node)[["label"]]
  node.graphic.attrs <- xmlAttrs(xmlElementsByTagName(el = current.node, name = "graphics")[[1]])
  nodes.x.pos[i] <- node.graphic.attrs[["x"]]
  nodes.y.pos[i] <- node.graphic.attrs[["y"]]
  nodes.color[i] <- node.graphic.attrs[["fill"]]
  
  node.attributes <- xmlElementsByTagName(el = current.node, name = "att")

  for(j in 1:length(attributes.names))
  {
    current.attribute <- xmlAttrs(node.attributes[[j]])
    attributes.df[i,current.attribute[["name"]]] <- current.attribute[["value"]]
  }
}

head(attributes.df)

nodes.df <- data.frame(names=nodes.names,x.pos=nodes.x.pos,y.pos=nodes.y.pos,color=nodes.color)
nodes.df <- cbind(nodes.df,attributes.df)
head(nodes.df)

## Remove selected column
colnames(nodes.df)
nodes.df <- nodes.df[,c(-5,-6,-7)]
head(nodes.df)

## Load network in gml format
library(igraph)
brc1.graph <- read.graph(file="BRC1_transcriptional_network.graphml", format = "graphml")
vertex.names <- as.vector(nodes.df$names)

## Store network connectivity
network.tfs <- read.table(file="network_tfs.txt",header=T,as.is=T)
network.adj <- as.data.frame(matrix(0,nrow=number.nodes,ncol=nrow(network.tfs)))
colnames(network.adj) <- network.tfs$AGI
is.data.frame(network.adj)


## Initialise vectors to store topological parameters
brc1.indegree <- vector(mode="numeric",length=length(vertex.names))
brc1.outdegree <- vector(mode="numeric",length=length(vertex.names))
brc1.trans <- vector(mode="numeric",length=length(vertex.names))
brc1.close <- vector(mode="numeric",length=length(vertex.names))
brc1.between <- vector(mode="numeric",length=length(vertex.names))
brc1.eccent <- vector(mode="numeric",length=length(vertex.names))

## Loop to retrieve regulators and topological parameters for each node.
for(i in 1:number.nodes)
{
  network.adj[i,neighbors(graph = brc1.graph, v=vertex.names[i], mode="in")$name] <- 1

  brc1.indegree[i] <- degree(graph = brc1.graph, v=vertex.names[i],mode = "in")
  brc1.outdegree[i] <- degree(graph = brc1.graph, v=vertex.names[i],mode = "out")
  brc1.trans[i] <- transitivity(graph = brc1.graph, type = "local", vids=vertex.names[i])
  brc1.close[i] <- closeness(graph = brc1.graph, vids=vertex.names[i], normalized = TRUE)
  brc1.between[i] <- betweenness(graph = brc1.graph, v = vertex.names[i], normalized = TRUE)
  brc1.eccent[i] <- eccentricity(graph = brc1.graph, v = vertex.names[i])
}

colnames(network.adj) <- network.tfs$name
head(network.adj)

topological.parameters <- data.frame(indegree = brc1.indegree,
                                     outdegree = brc1.outdegree,
                                     transitivity = brc1.trans,
                                     closeness = brc1.close,
                                     betweenness = brc1.between,
                                     eccentricity = brc1.eccent)
head(topological.parameters)
## Add connectivity info
nodes.df <- cbind(nodes.df,network.adj)
nodes.df <- cbind(nodes.df,topological.parameters)
head(nodes.df)

is.data.frame(nodes.df)

## Write network representation
write.table(nodes.df, 
            file="brc1_transcriptional_network.tsv", 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE)
network.data <- read.table(file="brc1_transcriptional_network.tsv",header = TRUE,as.is=TRUE,sep="\t",comment.char = "")
