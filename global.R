## Author:  Francisco J. Romero-Campero
## Contact: Francisco J. Romero-Campero - fran@us.es 

# Load neccesary libraries
library(shiny)
library(shinythemes)
library(ChIPpeakAnno)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(Biostrings)
library(seqinr)
library(org.At.tair.db)
library(igraph)
library(ggplot2)
library(ggrepel)
library(stringr)
library(clusterProfiler)
library(pathview)
library(shinycssloaders)
library(shinyWidgets)
library(shinyjs)
library(SuperExactTest)
library(VennDiagram)

##Load the network data
network.data <- read.table(file="data/brc1_network.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "",comment.char = "%",check.names=F)
rownames(network.data) <- network.data$names

## Tranforming coordinates for a better visualization
x.coord <- as.numeric(network.data$y.pos)
y.coord <- as.numeric(network.data$x.pos)

network.data$x.pos <- x.coord
network.data$y.pos <- y.coord

pos.data <- t(matrix(data=c(x.coord,y.coord),ncol=2))

## Transcription factors AGI ids and names
tfs.names <- c("ATAF1", "bZIP52", "PIF3", "MYB3", "ZAT10", 
               "ERF055", "VIP1", "ERF014", "NAC018", "NAP",
               "ANAC032", "RGA", "HB21", "HB6", "SVP",
               "ABI5", "IBH1", "BRC1", "ERF035", "GBF2",
               "PDF2", "HSFB2B", "TCX2", "WRKY18", "ABF3",
               "HB40", "ZAT6", "DREB2A", "MYB56", "ERF003",
               "NAC6/ORE1", "HAT2", "SPCH", "DOF5_4", "NAC102",
               "HB53", "HEC1","GBF3", "HSFB2A")

tf.ids <- c("AT1G01720", "AT1G06850" ,"AT1G09530", "AT1G22640", "AT1G27730",
            "AT1G36060", "AT1G43700", "AT1G44830", "AT1G52880", "AT1G69490",
            "AT1G77450", "AT2G01570", "AT2G18550", "AT2G22430", "AT2G22540",
            "AT2G36270", "AT2G43060", "AT3G18550", "AT3G60490", "AT4G01120",
            "AT4G04890", "AT4G11660", "AT4G14770", "AT4G31800", "AT4G34000",
            "AT4G36740", "AT5G04340", "AT5G05410", "AT5G17800", "AT5G25190",
            "AT5G39610", "AT5G47370", "AT5G53210", "AT5G60850", "AT5G63790",
            "AT5G66700", "AT5G67060", "AT2G46270", "AT5G62020")

names(tf.ids) <- tfs.names

tfs.order <- c("BRC1", "ABI5", "ABF3", "GBF3", "GBF2", "bZIP52",
               "ANAC032", "ATAF1", "NAC6/ORE1","NAC102", "NAP", "NAC018",
               "HB53","HB40","HB21","HAT2","HB6","PDF2",
               "MYB56", "MYB3",
               "HEC1", "SPCH", "IBH1", "PIF3",
               "ERF003", "DREB2A", "ERF035", "ERF014", "ERF055",
               "ZAT6", "ZAT10",
               "SVP", "WRKY18","HSFB2A","HSFB2B", "DOF5_4", "RGA", "TCX2","VIP1")

tf.ids <- tf.ids[tfs.order]
tfs.names <- tfs.order

tf.colors <- c("gold", "firebrick1", "firebrick2", "firebrick3", "firebrick", "firebrick4",
               "springgreen1","springgreen2","springgreen3","springgreen4","seagreen3","seagreen4",
               "salmon", "salmon1", "salmon2", "salmon3", "salmon4","brown",
               "red","red4",
               "deepskyblue1", "deepskyblue2", "deepskyblue3", "deepskyblue4",
               "plum", "plum1", "plum2", "plum3", "plum4",
               "peachpuff3", "peachpuff4",
               "darkgreen", "paleturquoise3", "rosybrown", "rosybrown4", "purple", "navyblue", "limegreen", "lightsteelblue3")
names(tf.colors) <- tfs.names

cluster.names <- c("Cluster UP_C1", "Cluster UP_C2", "Cluster UP_C3", "Cluster UP_C4",
                   "Cluster UP_C5", "Cluster UP_C6", "Cluster DOWN_C1", "Cluster DOWN_C2",
                   "Cluster DOWN_C3")

tfs.network.data <- subset(network.data, names %in% tf.ids) 

cluster.pos <- data.frame(x.pos = c(-1700,-900,-200,700,900,700,-100,-1100,-1800), 
                          y.pos = c(1200,1900,2050,1300,400,-1500,-2500,-2600,-1000),
                          cluster.name=cluster.names[1:9])

## Extract gene ids
genes <- sort(network.data$name)

## Load and extract Arabidopsis thaliana annotation regarding genes, exons and cds 
txdb <- TxDb.Athaliana.BioMart.plantsmart28
genes.data <- subset(genes(txdb), seqnames %in% c("1","2","3","4","5")) ## only nuclear genes are considered
genes.data <- as.data.frame(genes.data)
exons.data <- as.data.frame(exons(txdb))
cds.data <- as.data.frame(cds(txdb))

## Load all and circadian genes
alias2symbol.table <- AnnotationDbi::select(org.At.tair.db, 
                                            keys=keys(org.At.tair.db, keytype="ENTREZID"), 
                                            columns=c("SYMBOL", "TAIR"), keytype="ENTREZID")
ath.universe <- alias2symbol.table$TAIR
alias2symbol.table <- subset(alias2symbol.table, TAIR %in% genes)
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias[is.na(alias)] <- "" 
genes.selectize <- paste(names(alias), alias, sep=" - ")

## Setting conversion between alias and agis
agis <-alias2symbol.table$TAIR
names(agis) <- alias2symbol.table$SYMBOL
agis[is.na(agis)] <- ""

## Color vectors
line.colors <- c("blue","red", "darkgreen","black","#663300","#99003d","#b3b300","#4d0039","#4d2600","#006666","#000066","#003300","#333300","#660066")
area.colors <- c("skyblue","salmon", "lightgreen","lightgrey","#ffcc99","#ff99c2","#ffffb3","#ffe6f9","#ffe6cc","#80ffff","#b3b3ff","#99ff99","#e6e600","#ffb3ff")

## Load chromosome sequences
chr1 <- getSequence(read.fasta(file = "data/athaliana_genome/chr1.fa",seqtype = "AA"))[[1]] 
#chr1 <- getSequence(read.fasta(file = "data/athaliana_genome/chr1.fa",seqtype = "AA"))[[1]]
chr2 <- getSequence(read.fasta(file = "data/athaliana_genome/chr2.fa",seqtype = "AA"))[[1]]
chr3 <- getSequence(read.fasta(file = "data/athaliana_genome/chr3.fa",seqtype = "AA"))[[1]]
chr4 <- getSequence(read.fasta(file = "data/athaliana_genome/chr4.fa",seqtype = "AA"))[[1]]
chr5 <- getSequence(read.fasta(file = "data/athaliana_genome/chr5.fa",seqtype = "AA"))[[1]]

## Function to compute the reverse complement
reverse.complement <- function(dna.sequence)
{
  return(c2s(comp(rev(s2c(dna.sequence)),forceToLower = FALSE)))
}

## Load Position Weight Matrices
## Open file connection
con <- file("data/jaspar_motifs/pfm_plants_20180911.txt",open = "r")

## Empty list for storing PWM
motifs.pwm <- vector(mode="list",length = 453)
motif.ids <- vector(mode="character",length=453)
motif.names <- vector(mode="character",length=453)

## Load 64 PWM
for(j in 1:454)
{
  ## First line contains motif id and name
  first.line <- readLines(con,1)
  
  motif.ids[j] <- strsplit(first.line,split=" ")[[1]][1]
  motif.names[j] <- strsplit(first.line,split=" ")[[1]][2]
  
  ## Next four line contians probabilites for each nucleotide
  a.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  c.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  g.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  t.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  
  ## Construct PWM
  motif.pwm <- matrix(nrow = 4,ncol=length(a.row))
  
  motif.pwm[1,] <- a.row
  motif.pwm[2,] <- c.row
  motif.pwm[3,] <- g.row
  motif.pwm[4,] <- t.row
  
  rownames(motif.pwm) <- c("A","C","G","T")
  
  motifs.pwm[[j]] <- prop.table(motif.pwm,2)
}

## Close file connection
close(con)

## Naming list with PWM
names(motifs.pwm) <- motif.names
names(motif.ids) <- motif.names

## Load bed files for each transcription factor
bed.file.names <- c("data/bed_files/ABF3_ABA_optimal_narrowPeak_p16.bed",
                    "data/bed_files/ABI5_col_v3h.narrowPeak",
                    "data/bed_files/ANAC032_ABA_optimal_narrowPeak_p16.bed",
                    "data/bed_files/ANAC092_ORE1_col_v3a",
                    "data/bed_files/ANAC102_ABA_optimal_narrowPeak_p16.bed",
                    "data/bed_files/ATAF1.narrowPeak",
                    "data/bed_files/ATHB21 col_a.narrowPeak",
                    "data/bed_files/ATHB40_col_a.narrowPeak",
                    "data/bed_files/ATHB53_col_a.narrowPeak",
                    "data/bed_files/brc1_peaks.bed",
                    "data/bed_files/bZIP52.narrowPeak",
                    "data/bed_files/DREB2A_ABA_optimal_narrowPeak_p16.bed",
                    "data/bed_files/erf55.narrowPeak",
                    "data/bed_files/ERF14.narrowPeak",
                    "data/bed_files/ERF035 AT3G60490.narrowPeak",
                    "data/bed_files/ESE3.narrowPeak",
                    "data/bed_files/GBF2_ABA_optimal_narrowPeak_p16.bed",
                    "data/bed_files/HAT2.narrowPeak",
                    "data/bed_files/HB6_ABA_optimal_narrowPeak_p16.bed",
                    "data/bed_files/HSF6.narrowPeak",
                    "data/bed_files/HSF7.narrowPeak",
                    "data/bed_files/MYB3_ABA_optimal_narrowPeak_p16.bed",
                    "data/bed_files/MYB56.narrowPeak",
                    "data/bed_files/NAM.narrowPeak",
                    "data/bed_files/NAP.narrowPeak",
                    "data/bed_files/OBP4.narrowPeak",
                    "data/bed_files/PDF2.narrowPeak",
                    "data/bed_files/SVP.narrowPeak",
                    "data/bed_files/TCX2.narrowPeak",
                    "data/bed_files/VIP.narrowPeak",
                    "data/bed_files/WRKY18.narrowPeak",
                    "data/bed_files/ZAT6_ABA_optimal_narrowPeak_p16.bed",
                    "data/bed_files/ZAT10 STZ_col.narrowPeak",
                    "data/bed_files/hec1.narrowPeak",
                    "data/bed_files/PIF3.narrowPeak",
                    "data/bed_files/RGA.bed",
                    "data/bed_files/IBH1.bed",
                    "data/bed_files/SPCH.bed",
                    "data/bed_files/GBF3.bed")

## Load bed files
bed.files <- vector(mode = "list",length = length(bed.file.names))
for(i in 1:length(bed.file.names))
{
  bed.files[[i]] <- read.table(file=bed.file.names[i],header = F, as.is = T)
}

names(bed.files) <- c("ABF3", "ABI5", "ANAC032", "NAC6/ORE1", "NAC102", "ATAF1", "HB21", 
                      "HB40", "HB53", "BRC1", "bZIP52", "DREB2A", "ERF055", "ERF014", 
                      "ERF035", "ERF003", "GBF2", "HAT2", "HB6", "HSFB2A","HSFB2B", "MYB3", "MYB56", 
                      "NAC018", "NAP",  "DOF5_4", "PDF2", "SVP", "TCX2", "VIP1", "WRKY18", "ZAT6", 
                      "ZAT10", "HEC1", "PIF3","RGA","IBH1", "SPCH","GBF3")

## TF binding sites colors and symbol shapes
symbol.shapes <- c(17, 18, 19, 15)
symbol.color <- c("blue", "red", "darkgreen", "magenta")

## Node colors to represent in the global transcriptional network
node.colors <- network.data$color

## Extract adjacency matrix
adj.global.matrix <- as.matrix(network.data[,tfs.names])
rownames(adj.global.matrix) <- network.data$names

## Function to determine cluster membership
cluster.member <- function(clusters.vector,cluster)
{
  return(cluster %in% clusters.vector)
}

## Function to generate output table
create.output.table <- function(input.gene.df,alias,tfs.names)
{
  output.selected.genes.df <- data.frame(matrix(nrow=nrow(input.gene.df), ncol=6))
  colnames(output.selected.genes.df) <- c("AGI ID", "Gene Name", "Gene Description", "Regulators","log2FC","Cluster")
  output.selected.genes.df$`Gene Description` <- input.gene.df$annotation
  output.selected.genes.df$log2FC <- input.gene.df$FC
  output.selected.genes.df$Cluster <- input.gene.df$cluster
  i <- 1
  for(i in 1:nrow(output.selected.genes.df))
  {
    tair.link <- paste0("https://www.arabidopsis.org/servlets/TairObject?type=locus&name=",input.gene.df[i,1])
    output.selected.genes.df[i,1] <- paste(c("<a href=\"",
                                             tair.link,
                                             "\" target=\"_blank\">",
                                             input.gene.df[i,1], "</a>"),
                                           collapse="")
    output.selected.genes.df[i,2] <- alias[input.gene.df[i,1]]
    output.selected.genes.df[i,4] <- paste(tfs.names[which(input.gene.df[i,tfs.names] == 1)],collapse=", ")
  }
  
  return(output.selected.genes.df)
}

## Function to generate output table to download
create.downloadable.output.table <- function(input.gene.df,alias,tfs.names)
{
  output.selected.genes.df <- data.frame(matrix(nrow=nrow(input.gene.df), ncol=6))
  colnames(output.selected.genes.df) <- c("AGI ID", "Gene Name", "Gene Description", "Regulators","log2FC","Cluster")
  output.selected.genes.df$`Gene Description` <- input.gene.df$annotation
  output.selected.genes.df$log2FC <- input.gene.df$FC
  output.selected.genes.df$Cluster <- input.gene.df$cluster
  i <- 1
  for(i in 1:nrow(output.selected.genes.df))
  {
    output.selected.genes.df[i,1] <- input.gene.df[i,1]
    output.selected.genes.df[i,2] <- alias[input.gene.df[i,1]]
    output.selected.genes.df[i,4] <- paste(tfs.names[which(input.gene.df[i,tfs.names] == 1)],collapse=", ")
  }
  
  return(output.selected.genes.df)
}

## Auxiliary function to compute enrichments for GO table
compute.enrichments <- function(gene.ratios, bg.ratios)
{
  gene.ratios.eval <- sapply(parse(text=gene.ratios),FUN = eval)
  bg.ratios.eval <- sapply(parse(text=bg.ratios),FUN = eval)
  enrichments <- round(x=gene.ratios.eval/bg.ratios.eval,digits = 2)
  enrichments.text <- paste(enrichments, " (", gene.ratios, "; ", bg.ratios, ")",sep="")
  
  return(enrichments.text)  
}

#GO links and tair link functions
go.link <- function(go.term)
{
  link <- paste0("http://amigo.geneontology.org/amigo/term/", go.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           go.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

gene.link.function <- function(gene.name)
{
  tair.link <- paste(c("https://www.arabidopsis.org/servlets/TairObject?name=",
                       gene.name,
                       "&type=locus"),collapse="")
  gene.link <- paste(c("<a href=\"",
                       tair.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## KEGG pathway link
kegg.pathway.link <- function(kegg.pathway)
{
  link <- paste0("https://www.genome.jp/kegg-bin/show_pathway?",kegg.pathway)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           kegg.pathway, "</a>"),
                         collapse = "")
  return(complete.link)
}

## KEGG module link
kegg.module.link <- function(kegg.module)
{
  link <- paste0("https://www.genome.jp/kegg-bin/show_module?",kegg.module)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           kegg.module, "</a>"),
                         collapse = "")
  return(complete.link)
}

## TFBS jaspar link
tfbs.link <- function(motif.id)
{
  link <- paste0("http://jaspar.genereg.net/matrix/",motif.id)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           motif.id, "</a>"),
                         collapse = "")
  return(complete.link)
}
