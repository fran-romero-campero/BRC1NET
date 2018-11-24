# R script for pre-processing 
# Copyright (C) 2018  Francisco J. Romero-Campero, Pedro de los Reyes
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
# Date: August 2018

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

for(i in 1:number.nodes)
{
  current.node <- node.elements[[i]]
  nodes.names[i] <- xmlAttrs(current.node)[["label"]]
  node.graphic.attrs <- xmlAttrs(xmlElementsByTagName(el = current.node, name = "graphics")[[1]])
  nodes.x.pos[i] <- node.graphic.attrs[["x"]]
  nodes.y.pos[i] <- node.graphic.attrs[["y"]]
  
  node.attributes <- xmlElementsByTagName(el = current.node, name = "att")

  for(j in 1:length(attributes.names))
  {
    current.attribute <- xmlAttrs(node.attributes[[j]])
    attributes.df[i,current.attribute[["name"]]] <- current.attribute[["value"]]
  }
}

nodes.df <- data.frame(names=nodes.names,x.pos=nodes.x.pos,y.pos=nodes.y.pos)
head(nodes.df)

## Add info regarding clusters
peak.zt0 <- read.table(file="clusters/peak_ZT0.txt",as.is=T)[[1]]
peak.zt4 <- read.table(file="clusters/peak_ZT4.txt",as.is=T)[[1]]
peak.zt8 <- read.table(file="clusters/peak_ZT8.txt",as.is=T)[[1]]
peak.zt12 <- read.table(file="clusters/peak_ZT12.txt",as.is=T)[[1]]
peak.zt16 <- read.table(file="clusters/peak_ZT16.txt",as.is=T)[[1]]
peak.zt20 <- read.table(file="clusters/peak_ZT20.txt",as.is=T)[[1]]

cluster.genes <- c(peak.zt0, peak.zt4,peak.zt8,peak.zt12,peak.zt16,peak.zt20)
clusters.names <- paste0("peak", c(rep(0,length(peak.zt0)),
                                   rep(4,length(peak.zt4)),
                                   rep(8,length(peak.zt8)),
                                   rep(12,length(peak.zt12)),
                                   rep(16,length(peak.zt16)),
                                   rep(20,length(peak.zt20))))  

names(clusters.names) <- cluster.genes

cluster.classification <- clusters.names[as.vector(nodes.df$names)]
names(cluster.classification) <- NULL
nodes.df <- data.frame(nodes.df,cluster.classification)
write.table(x = nodes.df, file = "attractor_network.tsv",quote = FALSE,sep = "\t",row.names = FALSE)

