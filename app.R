## Author:  Francisco J. Romero-Campero
## Contact: Francisco J. Romero-Campero - fran@us.es 

## Input to test 
## input <- list(selected.multiple.tfs = c("BRC1", "HEC1", "ABI5"), cluster = "1")

## options(repos = BiocInstaller::biocinstallRepos())
## getOption("repos")

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
               "HB53", "HEC1","GBF3")

tf.ids <- c("AT1G01720", "AT1G06850" ,"AT1G09530", "AT1G22640", "AT1G27730",
            "AT1G36060", "AT1G43700", "AT1G44830", "AT1G52880", "AT1G69490",
            "AT1G77450", "AT2G01570", "AT2G02540", "AT2G22430", "AT2G22540",
            "AT2G36270", "AT2G43060", "AT3G18550", "AT3G60490", "AT4G01120",
            "AT4G04890", "AT4G11660", "AT4G14770", "AT4G31800", "AT4G34000",
            "AT4G36740", "AT5G04340", "AT5G05410", "AT5G17800", "AT5G25190",
            "AT5G39610", "AT5G47370", "AT5G53210", "AT5G60850", "AT5G63790",
            "AT5G66700", "AT5G67060", "AT2G46270")

names(tf.ids) <- tfs.names

tfs.order <- c("BRC1", "ABI5", "ABF3", "GBF3", "GBF2", "bZIP52",
               "ANAC032", "ATAF1", "NAC6/ORE1","NAC102", "NAP", "NAC018",
               "HB53","HB40","HB21","HAT2","HB6","PDF2",
               "MYB56", "MYB3",
               "HEC1", "SPCH", "IBH1", "PIF3",
               "ERF003", "DREB2A", "ERF035", "ERF014", "ERF055",
               "ZAT6", "ZAT10",
               "SVP", "WRKY18","HSFB2B", "DOF5_4", "RGA", "TCX2","VIP1")

tf.ids <- tf.ids[tfs.order]
tfs.names <- tfs.order

tf.colors <- c("gold", "firebrick1", "firebrick2", "firebrick3", "firebrick", "firebrick4",
               "springgreen1","springgreen2","springgreen3","springgreen4","seagreen3","seagreen4",
               "salmon", "salmon1", "salmon2", "salmon3", "salmon4","brown",
               "red","red4",
               "deepskyblue1", "deepskyblue2", "deepskyblue3", "deepskyblue4",
               "plum", "plum1", "plum2", "plum3", "plum4",
               "peachpuff3", "peachpuff4",
               "darkgreen", "paleturquoise3", "rosybrown", "purple", "navyblue", "limegreen", "lightsteelblue3")
names(tf.colors) <- tfs.names


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
                      "ERF035", "ERF003", "GBF2", "HAT2", "HB6", "HSFB2B", "MYB3", "MYB56", 
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


# Define UI for application that draws a histogram
ui <- fluidPage(
  ##shinythemes::themeSelector(),
  theme = shinytheme("sandstone"),
  ##theme = shinytheme("simplex"),
  ##theme = shinytheme("journal"),
  
  # Application title, introductory text and sidebar navigation
  fluidRow(
    column(
      width = 2,
      img(src='bud.jpg', align = "center", width=180),
      tags$br(),
      radioButtons(inputId = "navigation_bar", width="100%",selected="home",
                   label="",
                   choices=c(
                      "Home" = "home",
                     "Individual gene analysis" = "individual_gene",
                     "Multiple transcription factor analysis" = "multiple_gene",
                     "Tutorials" = "tutorials",
                     "GitHub repository" = "github",
                     "Citation and Contact" = "citation"
                   ))),
    column(
      width = 8,
      tags$div(align = "center", 
               tags$h1("The TCP transcription factor ",tags$i(tags$b("BRANCHED1,")),
                       " orchestrates several gene regulatory networks to promote axillary
                       bud dormancy in ", tags$i("Arabidopsis"))),
      tags$br(),tags$br(),
      conditionalPanel(condition = "input.navigation_bar == 'home'",
        tags$div(align="justify", "This online tool allows research to explore the transcriptional
          network downstream of", tags$i("BRC1"),". This network was constructed by combining genome-wide transcriptional 
          profiling of active and dormant buds, and of seedlings after BRC1 induction, with 
          determination of genome-wide BRC1 binding sites by Chromatin Immunoprecipitation 
          sequencing (ChIP-seq). We identified nine co-expressed gene clusters strongly dependent on 
          BRC1 function and, within these clusters, a group of genes encoding TFs which are 
          direct targets of BRC1. This specific set of regulators probably plays a key role 
          together with BRC1 in mediating amplification and maintenance of the observed 
          transcriptional responses triggered by BRC1. Use the navigation bar on the left to obtain
          insights into the molecular mechanisms that operate directly downstream of BRC1 to promote axillary 
          bud arrest.")
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'individual_gene'",
      tags$div(align="justify", tags$b("BRC1NET"), "allows researchers to explore the coordinated regulation of several
               BRC1 dependent transcription factors over an", tags$b("individually selected"), "gene. Follow the steps below:",
               tags$ol(
                 tags$li("Select a specific gene from our network using the", tags$b("Target Gene"),
                         "dropdown menu on the left below. You can enter either the AGI identifier or primary symbol for the",
                         "gene of interest."),
                 tags$li("Select several transcription factors (TFs) to explore their regulation over the
                         previously selected target gene using the",tags$b("Select Transcription Factors"), "checkboxes on the left below."),
                 tags$li("Results will be depicted below showing the selected target gene genomic location and the 
                          peaks detected in our analysis of the correspondig TFs ChIP-seq or DAP-seq data. 
                          Here you can specify the length of the gene promoter and 3' region as well as DNA 
                          TFs binding motifs to search for in the detected peak regions with the specified 
                          score for identity.")
               ))
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'multiple_gene'",
                       tags$div(align="justify", tags$b("BRC1NET"), "allows researchers to explore the coordinated regulation of several 
                                BRC1-dependent transcription factors or regulators over their common targets. Follow the steps below:", 
                                
                                tags$ol(
                                  tags$li("Select your TFs of interest using the", tags$b("Select Transcription Factors"),
                                          "checkbox menu on the left below."),
                                  tags$li("You can also select a gene cluster using the dropdown menu", tags$b("Select a gene cluster:")),
                                  tags$li("Check the box", tags$b("Visualize Edges"), "when you want to depict arrows from TFs to their target genes."),
                                  tags$li("Click on the", tags$b("SELECT GENES"), "to visualize your selected TFs common target genes located
                                           in the specified cluster. Explore the different tabs to 
                                           download a table with the selected genes, perform a signficance analysis of the overlap between the selected 
                                           TFs targets and the specified gene cluster as well as functional enrichment analysis.")
                                )
                       )
      )
      
    ),
    column(
      width = 2,
      img(src='cnb.jpg', align = "center", width=100),
      img(src='logo_ibvf.jpg', align = "center", width=100),
      img(src='logo_us.png', align = "center", width=100),
      tags$br(),tags$br(),tags$br(),
      img(src='logo_csic.jpg', align = "center", width=100)
    )
  ),
  
  ## Separation for different tools
  tags$br(),tags$br(),
  
  ## Conditional panel for individual gene analysis
  conditionalPanel(condition = "input.navigation_bar == 'individual_gene'",
                   fluidRow(
                     column(width = 3,
                            ## Select target gene to study
                            selectizeInput(inputId = "target.gene",
                                           label = "Target Gene:",
                                           choices = genes.selectize,
                                           multiple = FALSE),

                            ## Check box for the TF peaks to represent
                            checkboxGroupInput(inputId = "selected.tfs",
                                               label = "Select Transcription Factors:",
                                               choices = c("Select All Transcription Factors",
                                                           tfs.names),
                                               inline = TRUE,width = "100%")
                     ),

                     column(width = 9,
                                                 column(wellPanel(
                                                   ## Numeric input for promoter length
                                                   numericInput(inputId = "promoter.length",
                                                                label = "Promoter Length",
                                                                value = 2000,
                                                                min = 500,
                                                                max = 2000,
                                                                step = 100),
                                                   ## Numeric input for 3' length
                                                   numericInput(inputId = "threeprime.length",
                                                                label = "3' Length",
                                                                value = 500,
                                                                min = 100,
                                                                max = 500,
                                                                step = 100)),width=3),
                                                 column(wellPanel(
                                                   ## Selectize to choose target gene to represent
                                                   selectizeInput(inputId = "selected.motifs",
                                                                  label = "Select Motifs",
                                                                  selected =c("G-box","BRC1"),
                                                                  choices = motif.names,
                                                                  multiple = TRUE),
                                                   
                                                   ## Checkbox to select all available motifs
                                                   checkboxInput(inputId = "all.motifs",
                                                                 label = "Select All Motifs:",
                                                                 value = FALSE),

                                                   ## Numeric input for PWM score
                                                   numericInput(inputId = "min.score.pwm",
                                                                label = "Motif Identification Score:",
                                                                value = 100,
                                                                min = 80,
                                                                max = 100,
                                                                step = 5)),width=9),

                                                 fluidRow(
                                                   column(
                                                     plotOutput(outputId = "peak_plot"),
                                                     width=12)
                                                 )
                     )
                   )
  ),
  
  
  ## Conditional panel for multiple transcription factors and regulators analysis
  conditionalPanel(condition = "input.navigation_bar == 'multiple_gene'",
                   fluidRow(
                     column(width = 3,
                            ## Check box for the TFs to analyse
                            checkboxGroupInput(inputId = "selected.multiple.tfs",
                                               label = "Select Transcription Factors:",
                                               choices = tfs.names,
                                               inline = TRUE,width = "100%"),
                            tags$b("Select a gene cluster:"),
                            selectInput(inputId = "cluster", label="", 
                                        choices = c("No selected cluster" = "any", 
                                                    "Cluster 1" = "1", "Cluster 2" = "2", "Cluster 3" = "3",
                                                    "Cluster 4" = "4", "Cluster 5" = "5", "Cluster 6" = "6",
                                                    "Cluster 7" = "7", "Cluster 8" = "8", "Cluster 9" = "9"),
                                                  selected = "any", multiple = FALSE, selectize = TRUE),
                            checkboxInput(inputId =  "edges",label = "Visualize Edges",value = FALSE),
                            
                            radioButtons(inputId = "selection.mode",label = "Mode of selection:",
                                         choices = c("Common gene targets","All gene targets"),selected = "Common gene targets",inline = F),
                            
                            actionButton(inputId = "go_multiple",label = "Select Genes")
                     ),
                     
                     column(width = 9,
                            tabsetPanel(type = "tabs",
                                        tabPanel(title = "Network Visualization",
                                                 plotOutput("networkPlot")
                                        ),
                                        tabPanel(title = "Gene Table",
                                                 dataTableOutput(outputId = "outputTable"),
                                                 uiOutput(outputId = "download_ui_for_table")
                                        ),
                                        tabPanel(title = "Overlap Significance",
                                                 tags$br(),
                                                 tags$div(align="justify", "In this section, we present the results of a significance analysis of the
                                                           overlap between the targets of the selected transcription factors and gene
                                                           with a specific expresion pattern."),
                                                 tags$br(),
                                                 textOutput("overlap.message"),
                                                 textOutput("overlap.significance.text"),
                                                 tags$br(),
                                                 tags$br(),
                                                 tags$div(align="center",
                                                          plotOutput("venn.diagram.plot"))
                                        ),
                                        tabPanel(title = "Functional Enrichment",
                                                 tabsetPanel(type = "tabs",
                                                             tabPanel(title = "GO Enrichment",
                                                                      tags$br(),
                                                                      tags$div(align="justify", "In this section you can perform a GO term
                                                                               enrichment analysis over the selected genes. First
                                                                               of all, you need to choose the background set of genes between
                                                                               the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in ATTRACTOR:"),
                                                                      tags$br(),
                                                                      radioButtons(inputId = "go.background", width="100%",selected="allgenome",
                                                                                   label="",
                                                                                   choices=c(
                                                                                     "Complete genome" = "allgenome",
                                                                                     "Genes in network" = "onlynet"
                                                                                   )), tags$br(),
                                                                      actionButton(inputId = "goterm",label = "GO terms analysis"),
                                                                      tags$br(),
                                                                      tags$br(),
                                                                      shinyjs::useShinyjs(),
                                                                      hidden(div(id='loading.div',h3('Please be patient, computing GO enrichment ...'))),
                                                                      tags$br(),
                                                                      tabsetPanel(type = "tabs",
                                                                                  tabPanel(title = "GO table",
                                                                                           tags$br(), tags$br(),
                                                                                           htmlOutput(outputId = "textGOTable"),
                                                                                           tags$br(), tags$br(),
                                                                                           dataTableOutput(outputId = "output_go_table"),
                                                                                           htmlOutput(outputId = "revigo"),
                                                                                           uiOutput(outputId = "download_ui_for_go_table")
                                                                                  ),
                                                                                  tabPanel(title = "GO map",
                                                                                           tags$br(), tags$br(),
                                                                                           htmlOutput(outputId = "gomap_text"),
                                                                                           tags$br(),
                                                                                           div(style= "overflow:scroll; height:500px;",
                                                                                               addSpinner(plotOutput("gomap"), spin = "circle", color = "#E41A1C"))),
                                                                                  tabPanel(title = "GO barplot",
                                                                                           tags$br(),
                                                                                           htmlOutput(outputId = "barplot_text"),
                                                                                           tags$br(),
                                                                                           addSpinner(plotOutput("bar.plot"), spin = "circle", color = "#E41A1C")),
                                                                                  tabPanel(title = "GO Enrichment Map",
                                                                                           tags$br(), 
                                                                                           htmlOutput(outputId = "emapplot_text"),
                                                                                           tags$br(),
                                                                                           addSpinner(plotOutput(outputId = "emap.plot",inline=TRUE))),
                                                                                  tabPanel(title = "GO concept network",
                                                                                           tags$br(), 
                                                                                           htmlOutput(outputId = "cnetplot_text"),
                                                                                           tags$br(),
                                                                                           addSpinner(plotOutput(outputId = "cnet.plot",inline=TRUE)))
                                                                      )
                                                                      ),
                                                             tabPanel(title = "KEGG Pathway Enrichment",
                                                                      tags$br(),
                                                                      tags$div(align="justify", "In this section you can perform a KEGG pathways and modules
                                                                               enrichment analysis over the selected genes. First
                                                                               of all, you need to choose the background set of genes between
                                                                               the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in ATTRACTOR:"),
                                                                      tags$br(),
                                                                      radioButtons(inputId = "pathway_background", width="100%",selected="allgenome",
                                                                                   label="",
                                                                                   choices=c(
                                                                                     "Complete genome" = "allgenome",
                                                                                     "Genes in network" = "onlynet"
                                                                                   )),
                                                                      actionButton(inputId = "pathway_button",label = "KEGG pathway analysis"),
                                                                      tags$br(),
                                                                      tags$br(),
                                                                      shinyjs::useShinyjs(),
                                                                      hidden(div(id='loading.div.kegg',h3('Please be patient, computing KEGG pathway enrichment ...'))),
                                                                      tags$br(),
                                                                      tabsetPanel(type = "tabs",
                                                                                  tabPanel(title = "Enriched Pathway Table",
                                                                                           tags$br(), tags$br(),
                                                                                           htmlOutput(outputId = "no_kegg_enrichment"),
                                                                                           dataTableOutput(outputId = "output_pathway_table"),
                                                                                           uiOutput(outputId = "download_ui_for_kegg_table")
                                                                                  ),
                                                                                  tabPanel(title = "Enriched Pathway Visualization",
                                                                                           uiOutput(outputId = "kegg_selectize"),
                                                                                           imageOutput("kegg_image"),
                                                                                           br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                                                                           br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                                                                           br(), br(), br(), br(), br()
                                                                                  ),
                                                                                  tabPanel(title = "Enriched Module Table",
                                                                                           htmlOutput(outputId = "text_module_kegg"),
                                                                                           br(), br(),
                                                                                           dataTableOutput(outputId = "output_module_table")
                                                                                  )
                                                                      )
                                                                      )
                                                 )
                                        )
                             )
                     )
                                                 
                   )
  )

)

## ATTRACTOR server
server <- function(input, output, session) {
  
  # ## Select all code
  # observe({
  #   if ("Select All Transcription Factors" %in% input$selected.tfs) {
  #     # choose all the choices _except_ "Select All"
  #     updateSelectInput(session, "selected.tfs", selected = tfs.names)
  #   }
  # })

  ## Peak visualizer code
  output$peak_plot <- renderPlot({

    ## Select all code
    if ("Select All Transcription Factors" %in% input$selected.tfs) 
    {
      # choose all the choices _except_ "Select All"
      updateSelectInput(session, "selected.tfs", selected = tfs.names)
      currently.selected.tfs <- tfs.names
    } else
    {
      currently.selected.tfs <- input$selected.tfs
    }
     
      ## Sanity checks
      validate(
        need(length(input$selected.tfs) > 0 , "Please select a set of transcription factors"),
        need(input$target.gene, "Please select a target gene")
      )
      
      ## Extract target gene annotation 
      gene.name <-  strsplit(input$target.gene,split=" - ")[[1]][1]
      
      target.gene.body <- genes.data[gene.name,]
      target.gene.chr <- as.character(target.gene.body$seqnames)
      target.gene.start <- target.gene.body$start
      target.gene.end <- target.gene.body$end
      
      target.gene.strand <- as.character(target.gene.body$strand)
      
      ## Extract cds annotation
      cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      ## Extract exons annotation
      exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      ## Determine the genome range to plot including promoter, gene body and 3' UTR
      ## This depends on whether the gene is on the forward or reverse strand
      range.to.plot <- target.gene.body
      
      if(target.gene.strand == "+")
      {
        range.to.plot$start <- range.to.plot$start - input$promoter.length
        range.to.plot$end <- range.to.plot$end + input$threeprime.length
      } else if (target.gene.strand == "-")
      {
        range.to.plot$end <- range.to.plot$end + input$promoter.length
        range.to.plot$start <- range.to.plot$start - input$threeprime.length
      }
      
      ## Compute the length of the genome range to represent
      current.length <- range.to.plot$end - range.to.plot$start
      
      ## Determine upper limit of the graph
      number.tfs <- length(currently.selected.tfs) #length(input$selected.tfs)
      upper.lim <- 25 * length(currently.selected.tfs) #length(input$selected.tfs)
      
      ## Draw DNA strand
      gene.height <- -25
      cord.x <- 1:current.length
      
      plot(cord.x, rep(gene.height,length(cord.x)),type="l",col="black",lwd=3,ylab="",
           cex.lab=2,axes=FALSE,xlab="",main="",cex.main=2,
           ylim=c(-30,25 * length(currently.selected.tfs)), #length(input$selected.tfs)
           xlim=c(-3000,max(cord.x)))
      
      ## Extract exons for target gene
      exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      ## Transform exon coordinates to current range
      min.pos <- min(exons.data.target.gene$start)
      
      if(target.gene.strand == "+")
      {
        exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$promoter.length
        exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$promoter.length
      } else if(target.gene.strand == "-")
      {
        exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$threeprime.length
        exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$threeprime.length
      }
      
      ## Represent exons
      exon.width <- 2
      for(i in 1:nrow(exons.data.target.gene))
      {
        # Determine start/end for each exon
        current.exon.start <- exons.data.target.gene$start[i]
        current.exon.end <- exons.data.target.gene$end[i]
        
        ## Determine coordinates for each exon polygon and represent it
        exon.x <- c(current.exon.start,current.exon.end,current.exon.end,current.exon.start)
        exon.y <- c(gene.height + exon.width, gene.height + exon.width, gene.height - exon.width, gene.height - exon.width)
        
        polygon(x = exon.x, y = exon.y, col = "blue",border = "blue")
      }
      
      ## Extract cds for target gene
      cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      ## Transform cds coordinates to current range
      if(target.gene.strand == "+")
      {
        cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$promoter.length
        cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$promoter.length
      } else if (target.gene.strand == "-")
      {
        cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$threeprime.length
        cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$threeprime.length
      }
      
      cds.width <- 3
      for(i in 1:nrow(cds.data.target.gene))
      {
        # Determine current cds start/end
        current.cds.start <- cds.data.target.gene$start[i]
        current.cds.end <- cds.data.target.gene$end[i]
        
        # Determine curret cds coordinates for the polygon and represent it
        cds.x <- c(current.cds.start,current.cds.end,current.cds.end,current.cds.start)
        cds.y <- c(gene.height + cds.width, gene.height + cds.width, gene.height - cds.width, gene.height - cds.width)
        
        polygon(x = cds.x, y = cds.y, col = "blue",border = "blue")
      }
      
      ## Draw arrow to represent transcription direction 
      if(target.gene.strand == "+")
      {
        lines(c(input$promoter.length,input$promoter.length,input$promoter.length+100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
        lines(c(input$promoter.length+50,input$promoter.length+100),y=c(gene.height+6,gene.height+5),lwd=3)
        lines(c(input$promoter.length+50,input$promoter.length+100),y=c(gene.height+4,gene.height+5),lwd=3)
      } else if (target.gene.strand == "-")
      {
        lines(c(current.length - input$promoter.length, current.length - input$promoter.length, current.length - input$promoter.length-100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
        lines(c(current.length - input$promoter.length-50, current.length - input$promoter.length - 100),y=c(gene.height + 6, gene.height + 5),lwd=3)
        lines(c(current.length - input$promoter.length-50, current.length - input$promoter.length - 100),y=c(gene.height + 4, gene.height + 5),lwd=3)
      }
      
      ## Draw promoter range
      if(target.gene.strand == "+")
      {
        axis(side = 1,labels = c(- input$promoter.length, - input$promoter.length / 2,"TSS"),at = c(1,input$promoter.length/2,input$promoter.length),lwd=2,cex=1.5,las=2,cex=2)
      } else if(target.gene.strand == "-")
      {
        axis(side = 1,labels = c("TSS",- input$promoter.length / 2,- input$promoter.length),at = c(current.length-input$promoter.length,current.length-input$promoter.length/2, current.length),lwd=2,cex=1.5,las=2,cex=2)
      }
      
      selected.bed.files <- bed.files[currently.selected.tfs]#input$selected.tfs]
      
      ## Since ChIPpeakAnno needs more than one region to plot our region
      ## is duplicated 
      regions.plot <- GRanges(rbind(range.to.plot,range.to.plot))
      

      ## Draw peak regions for each TF and determing TF binding sequences
      
      ## Determine TFBS motifs to search for and Selecting an
      ## example motif if the user does not select any of them
      if(input$all.motifs)
      {
        selected.motifs.pwm <- motifs.pwm
      } else if(length(input$selected.motifs)==0)
      {
        selected.motifs.pwm <- motifs.pwm[c("G-box","BRC1")]
      }else
      {
        selected.motifs.pwm <- motifs.pwm[input$selected.motifs]
      }
      
      selected.motif.names <- names(selected.motifs.pwm)
      selected.motif.ids <- motif.ids[selected.motif.names]
      
      ## Initialize data frame containing TF binding sequences in the peak regions
      df.hits <- data.frame(0,0,"","","")
      colnames(df.hits) <- c("tf_number","position","id","name","seq")
      
      ## Width of the rectangle representing the peak region
      peak.width <- 1
      for(i in 1:length(currently.selected.tfs))#length(input$selected.tfs))
      {
        current.tf <- currently.selected.tfs[i] #input$selected.tfs[i]
        ## Extract bed file name 1 and read it
        current.peaks <- bed.files[[current.tf]]
        #current.peaks <- read.table(file=current.bed.file,header = F, as.is = T)
        peak.coordinates <- subset(current.peaks, V1 == range.to.plot$seqnames & V2 >= range.to.plot$start & V3 <= range.to.plot$end) 
        current.peaks.to.plot <- peak.coordinates[,2:3]
        
        ## Transform coordinates 
        current.peaks.to.plot <- current.peaks.to.plot - range.to.plot$start
        
        ## Check if there are peaks for the target gene
        if(nrow(current.peaks.to.plot) > 0)
        {
          ## Normalization
          #chip.signal.means[i, ] <- 10 * chip.signal.means[i, ] / max(chip.signal.means[i, ])
          
          #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
          for(j in 1:nrow(current.peaks.to.plot))
          {
            ## Extract start and end point of each peak region
            current.peak.start <- current.peaks.to.plot[j,1]
            current.peak.end <- current.peaks.to.plot[j,2]
            
            ## Computer coordinates for polygon and draw it
            peak.x <- c(current.peak.start,current.peak.end,
                        current.peak.end,current.peak.start)
            peak.y <- c(25*(i - 1) - 5 + peak.width, 25*(i - 1) - 5 + peak.width, 
                        25*(i - 1) - 5 - peak.width, 25*(i - 1) - 5 - peak.width)  
            
            #polygon(x = peak.x, y = peak.y, col = area.colors[i], border = line.colors[i],lwd=2)
            polygon(x = peak.x, y = peak.y, col = tf.colors[current.tf], border = tf.colors[current.tf],lwd=2)
                        
            ## Identify TF binding DNA motifs 
            peak.chr <- peak.coordinates[j, 1]
            peak.start <- peak.coordinates[j, 2]
            peak.end <- peak.coordinates[j, 3]
            
            ## Extract peak sequence
            if(peak.chr == "1")
            {
              peak.sequence <- c2s(chr1[peak.start:peak.end])
            } else if(peak.chr == "2")
            {
              peak.sequence <- c2s(chr2[peak.start:peak.end])
            } else if(peak.chr == "3")
            {
              peak.sequence <- c2s(chr3[peak.start:peak.end])
            } else if(peak.chr == "4")
            {
              peak.sequence <- c2s(chr4[peak.start:peak.end])
            } else if(peak.chr == "5")
            {
              peak.sequence <- c2s(chr5[peak.start:peak.end])
            }
            
            peak.rev.comp.sequence <- reverse.complement(peak.sequence)
            
            for(k in 1:length(selected.motifs.pwm))
            {
              motif.pwm <- selected.motifs.pwm[[k]]
              
              hits.fw <- matchPWM(motif.pwm, peak.sequence, 
                                  min.score = paste0(input$min.score.pwm,"%"))
              hits.fw.seqs <- as.data.frame(hits.fw)[[1]]
              hits.fw <- as(hits.fw, "IRanges")
              hits.fw.start <- start(hits.fw)
              hits.fw.end <- end(hits.fw)
              
              if(length(hits.fw.start) > 0)
              {
                df.hits.fw <- data.frame(rep(i,length(hits.fw.start)),
                                         ((hits.fw.start+hits.fw.end)/2) + current.peak.start,
                                         rep(selected.motif.ids[k],length(hits.fw.start)),
                                         rep(selected.motif.names[k],length(hits.fw.start)),
                                         hits.fw.seqs)
                colnames(df.hits.fw)  <- c("tf_number","position","id","name","seq")
                df.hits <- rbind(df.hits,df.hits.fw)
              }
              
              hits.rev <- matchPWM(motif.pwm, peak.rev.comp.sequence, 
                                   min.score = paste0(input$min.score.pwm,"%"))
              hits.rev.seqs <- as.data.frame(hits.rev)[[1]]
              hits.rev.seqs <- sapply(hits.rev.seqs,reverse.complement)
              names(hits.rev.seqs) <- NULL
              
              hits.rev <- as(hits.rev, "IRanges")
              hits.rev.start <- nchar(peak.sequence) - end(hits.rev) + 1
              hits.rev.end <- nchar(peak.sequence) - start(hits.rev) + 1
              
              if(length(hits.rev.start) > 0)
              {
                df.hits.rev <- data.frame(rep(i,length(hits.rev.start)),
                                          ((hits.rev.start+hits.rev.end)/2) + current.peak.start,
                                          rep(selected.motif.ids[k],length(hits.rev.start)),
                                          rep(selected.motif.names[k],length(hits.rev.start)),
                                          hits.rev.seqs)
                colnames(df.hits.rev)  <- c("tf_number","position","id","name","seq")
                df.hits <- rbind(df.hits,df.hits.rev)
              }
              
            }
            
          }
        }
      }
      
      ## Remove first line of the data frame added just for technical reason
      df.hits <- df.hits[-1,]
      
      ## Draw TF binding sites
      detected.tfbs <- unique(as.vector(df.hits$name))
      
      number.of.shapes <- ceiling(length(detected.tfbs) / length(symbol.color))
      
      necessary.shapes <- rep(symbol.shapes[1:number.of.shapes],each = length(detected.tfbs)/number.of.shapes)
      necessary.colors <- rep(symbol.color,number.of.shapes)
      
      if(length(detected.tfbs) > 0)
      {
        for(i in 1:length(detected.tfbs))
        {
          current.tfbs <- detected.tfbs[i]
          current.shape <- necessary.shapes[i]
          current.color <- necessary.colors[i]
          
          positions <- subset(df.hits, name == current.tfbs)
          
          for(j in 1:nrow(positions))
          {
            tf.to.draw <- positions$tf_number[j]
            pos.to.draw <- positions$position[j]
            
            points(x = pos.to.draw, y = 25*(tf.to.draw - 1) - 5 - 5*peak.width,
                   pch = current.shape, col = current.color, cex = 1)
          }
        }
        
        ## Add legend for TFBS
        legend.step <- 10
        for(i in 1:length(detected.tfbs))
        {
          points(x = -3000, y = upper.lim - (i-1)*legend.step, 
                 pch=necessary.shapes[i], col = necessary.colors[i],cex = 1)
          
          
          current.seq <- as.character(subset(df.hits,name == detected.tfbs[i])[["seq"]][[1]])
          current.label <- paste(c(detected.tfbs[i], "  -  ", current.seq ),collapse="")
          
          text(x = -2900, y = upper.lim - (i-1)*legend.step, labels = current.label,
               adj = 0,cex = 0.7)
        }
      }
      
      ## Draw profiles for TF binding
      for(i in 1:length(currently.selected.tfs))#length(input$selected.tfs))
      {
        #text(x = -50,y = 25*(i-1),labels = input$selected.tfs[i],adj = 1,col = tf.colors[input$selected.tfs[i]],font = 2)
        text(x = -50,y = 25*(i-1),labels = currently.selected.tfs[i],adj = 1,col = tf.colors[currently.selected.tfs[i]],font = 2)
      }
      
    })

  ## Multiple transcription factor code
  
  ## Initial/default visualization of ATTRACTOR
  output$networkPlot <- renderPlot({
    ggplot(network.data, aes(x.pos,y.pos)) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank()) + 
      geom_point(color=node.colors,size=1)
  },height = 700)
  
  ## Determine common targets and perform analysis when button is clicked
  observeEvent(input$go_multiple, {

    ## Determine targets of selected TFs
    selected.tfs.agi <- tf.ids[input$selected.multiple.tfs]
    
    #selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    #selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,input$selected.multiple.tfs] != 0)
    
    if(length(selected.tfs.agi) > 1)
    {
      if(input$selection.mode == "Common gene targets")
      {
        gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.agi)
      } else
      {
        gene.selection <- rowSums(selected.tfs.adj) >= 1
      }
      
    } else if (length(selected.tfs.agi) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }

    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]    
    
    selected.genes.df$cluster
    
    if(input$cluster != "any")
    {
      selected.genes.df <- subset(selected.genes.df, 
                                  sapply(strsplit(x = selected.genes.df$cluster,split = "|"), cluster.member,cluster = input$cluster))
    } 
    
    ## Node colors for representation
    selected.nodes.colors <- selected.genes.df$color
    
    ## Determine significance of overlap
    ## Number of sets to overlap
    if(input$cluster != "any")
    {
      number.of.sets <- length(input$selected.multiple.tfs) + 1
    } else
    {
      number.of.sets <- length(input$selected.multiple.tfs)
    }
    
    list.of.sets <- vector(mode = "list",length=number.of.sets)
    
    ## TFs targets
    for(i in 1:length(selected.tfs.agi))
    {
      list.of.sets[[i]] <- rownames(network.data)[which(network.data[,input$selected.multiple.tfs[[i]]] != 0)]
    }
    
    ## Gene cluster
    if(input$cluster != "any")
    {
      list.of.sets[[length(list.of.sets)]] <- rownames(subset(network.data, 
                                                       sapply(strsplit(x = network.data$cluster,split = "|"), cluster.member,cluster = input$cluster)))
      
      name.of.sets <- c(input$selected.multiple.tfs, input$cluster)
    } else
    {
      name.of.sets <- input$selected.multiple.tfs
    }
    
    names(list.of.sets) <- name.of.sets
    
    if(number.of.sets > 1)
    {
      ## Compute overlap p value and enrichment
      overlap.output <- supertest(x = list.of.sets, n = nrow(network.data))
      overlap.results <- summary(overlap.output)
      overlap.p.value <- tail(overlap.results$P.value, n=1)
      overlap.p.value <- round(overlap.p.value,digits = -log10(overlap.p.value)+2)
      overlap.enrichment <- round((overlap.results$Table)[["FE"]][nrow(overlap.results$Table)],digits=2)
      
      ## Ouput text for the overlap
      overlap.message <- "The overlap between the targets of the transcription factors"
      text.first.tfs <- paste(input$selected.multiple.tfs[1:(length(input$selected.multiple.tfs)-1)],collapse=",")
      text.all.tfs <- paste(c(text.first.tfs, "and", input$selected.multiple.tfs[length(input$selected.multiple.tfs)]),collapse=" ")
      overlap.message <- paste(overlap.message,text.all.tfs,sep=" ")
    
      if(input$cluster != "any")
      {
        overlap.message <- paste(overlap.message, paste(c("and the genes in cluster", 
                                                        input$cluster), collapse=" "),sep= " ")
      } 
    
      if(overlap.p.value < 0.05)
      {
        overlap.message <- paste(overlap.message, paste0(paste(c("is significant with a p-value of",overlap.p.value,"and an enrichment of", overlap.enrichment),collapse=" "),"."), sep=" ")
      } else
      {
        overlap.message <- paste(overlap.message, paste0(paste(c("is NOT significant according to a p-value of",overlap.p.value),collapse=" "),"."), sep = " ")
      }
      
    
      output$overlap.significance.text <- renderText(expr = {
        overlap.message
      })
    }
    
    ## Draw venn diagram
    if(number.of.sets == 1)
    {
      output$overlap.significance.text <- renderText(expr = {
        "Only a single transcription factor was selected. The analysis of 
        overlap significance is not applicable."
      })
    }
    if(number.of.sets == 2)
    {
      output$venn.diagram.plot <- renderPlot(expr = {
        grid.newpage()
        draw.pairwise.venn(area1 = length(list.of.sets[[1]]),
                           area2 = length(list.of.sets[[2]]),
                           cross.area = length(intersect(list.of.sets[[1]],list.of.sets[[2]])),
                           category = name.of.sets,
                           scaled = TRUE,lwd = 2,col = "black",
                           fill = c("blue","red"),alpha = 0.7,
                           cex = 2,cat.cex = 1.5,cat.pos = c(0,0))
      },width = 600,height = 600)
    } else if (number.of.sets == 3)
    {
      output$venn.diagram.plot <- renderPlot(expr = {
        grid.newpage()
        draw.triple.venn(area1 = length(list.of.sets[[1]]),
                         area2 = length(list.of.sets[[2]]),
                         area3 = length(list.of.sets[[3]]),
                         n12 = length(intersect(list.of.sets[[1]],list.of.sets[[2]])),
                         n23 = length(intersect(list.of.sets[[2]],list.of.sets[[3]])),
                         n13 = length(intersect(list.of.sets[[1]],list.of.sets[[3]])), 
                         n123 =  length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]),list.of.sets[[3]])),
                         category = name.of.sets,
                         scaled = TRUE,lwd = 2,col = "black",
                         fill = c("blue","red","darkgreen"),alpha = 0.7,
                         cex = 1.5,cat.cex = 1.5,cat.pos = c(-30,30,180))
      }, width = 600, height = 600)
    } else if (number.of.sets == 4)
    {
      output$venn.diagram.plot <- renderPlot(expr = {
        grid.newpage()
        draw.quad.venn(area1 = length(list.of.sets[[1]]),
                       area2 = length(list.of.sets[[2]]),
                       area3 = length(list.of.sets[[3]]),
                       area4 = length(list.of.sets[[4]]),
                       n12 = length(intersect(list.of.sets[[1]],list.of.sets[[2]])),
                       n13 = length(intersect(list.of.sets[[1]],list.of.sets[[3]])), 
                       n14 = length(intersect(list.of.sets[[1]],list.of.sets[[4]])),
                       n23 = length(intersect(list.of.sets[[2]],list.of.sets[[3]])),
                       n24 = length(intersect(list.of.sets[[2]],list.of.sets[[4]])),
                       n34 = length(intersect(list.of.sets[[3]],list.of.sets[[4]])),
                       n123 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]),list.of.sets[[3]])), 
                       n124 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]),list.of.sets[[4]])),
                       n134 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[3]]),list.of.sets[[4]])),
                       n234 = length(intersect(intersect(list.of.sets[[2]],list.of.sets[[3]]),list.of.sets[[4]])),
                       n1234 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]), intersect(list.of.sets[[3]],list.of.sets[[4]]))),
                       category = name.of.sets,
                       scaled = TRUE,lwd = 2,col = "black",
                       fill = c("blue","red","darkgreen","orange"),alpha = 0.7,
                       cex = 1.5,cat.cex = 1.5,cat.pos = c(0,0,0,0))
      }, width = 600, height = 600)
    } else if (number.of.sets > 4)
    {
      output$venn.diagram.plot <- renderPlot(expr = {
            plot(overlap.output, Layout = "landscape")
      })
    }
    
    ## Target gene representation on the network
    network.representation <- ggplot(network.data, aes(x.pos,y.pos)) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank()) + 
      geom_point(color=node.colors,size=1) +
      geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    
    ## Add edges to the network when selected
    if(input$edges)
    {
      for(i in 1:length(input$selected.multiple.tfs))
      {
        tf.xpos <- subset(network.data, names == tf.ids[input$selected.multiple.tfs[i]])[["x.pos"]]
        tf.ypos <- subset(network.data, names == tf.ids[input$selected.multiple.tfs[i]])[["y.pos"]]
        network.representation <- network.representation +
          annotate("segment",
                   x=rep(tf.xpos,nrow(selected.genes.df)),
                   y=rep(tf.ypos,nrow(selected.genes.df)),
                   xend=selected.genes.df$x.pos,
                   yend=selected.genes.df$y.pos, 
                   color="grey", arrow=arrow(type="closed",length=unit(0.1, "cm")))
      }
      
      network.representation <- network.representation + 
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    }
    
    ## Update network representation on the app
    output$networkPlot <- renderPlot({
      network.representation
    },height = 700)
    
    ## Output table with gene info
    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)
    
    ## Generate UI to download table and creating a downlodable table
    output$download_ui_for_table<- renderUI(
      tagList(downloadButton(outputId= "downloadData", "Download Selected Genes"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
    )
    
    output$downloadData<- downloadHandler(
      filename= function() {
        paste0(paste(input$selected.tfs,collapse = "_"), ".tsv")
      },
      content= function(file) {
        write.table(create.downloadable.output.table(input.gene.df=selected.genes.df,alias,tfs.names), 
                    file=file, 
                    sep = "\t", 
                    quote = FALSE,
                    row.names = FALSE)
      })
    
    
    
  })
  
  ##Perform GO terms enrichment analysis when button is clicked
  observeEvent(input$goterm,{
    
    ## Sanity checks
    validate(
      need(input$selected.multiple.tfs, "Please select a set of transcription factors")
    )

    ## Determine targets of selected TFs
    # selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    # selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,input$selected.multiple.tfs] != 0)
    
    if(length(input$selected.multiple.tfs) > 1)
    {
      
      if(input$selection.mode == "Common gene targets")
      {
        gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.agi)
      } else
      {
        gene.selection <- rowSums(selected.tfs.adj) >= 1
      }
      
    } else if (length(input$selected.multiple.tfs) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]    
    
    ## GO enrichment analysis
    
    ## Set the background to perform the GO terms enrichment analysis depending on the user selection
    if (input$go.background == "allgenome")
    {
      go.universe <- ath.universe
    } else
    {
      go.universe <- network.data$name
    }
    
    
    ## Show element when GO term enrichment analysis starts
    shinyjs::showElement(id = 'loading.div')

    enrich.go <- enrichGO(gene          = selected.genes.df$name,
                          universe      = go.universe,
                          OrgDb         = org.At.tair.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = FALSE,
                          keyType = "TAIR")

    ## Hide loading element when GO term enrichment is done
    shinyjs::hideElement(id = 'loading.div')
    
    ## Generate ouput table
    enrich.go.result <- as.data.frame(enrich.go)
    
    if(nrow(enrich.go.result) > 0)
    {
      ## GO term Description P-value Q-value Enrichment (SetRatio, BgRatio) Genes
      go.term.enrichments <- compute.enrichments(gene.ratios = enrich.go.result$GeneRatio,
                                                 bg.ratios = enrich.go.result$BgRatio)
      
      go.result.table <- data.frame(enrich.go.result$ID, enrich.go.result$Description,
                                    enrich.go.result$pvalue, enrich.go.result$qvalue,
                                    go.term.enrichments, 
                                    gsub(pattern = "/",replacement = " ",x = enrich.go.result$geneID),
                                    stringsAsFactors = FALSE)
      
      colnames(go.result.table) <- c("GO ID", "Description", "p-value", "q-value",
                                     "Enrichment (Target Ratio; BG Ration)","Genes")
      
      go.result.table.with.links <- go.result.table
      ## Add links to the genes
      genes.in.go.enrichment <- go.result.table$Genes
      
      ## Add link to genes
      for(i in 1:length(genes.in.go.enrichment))
      {
        go.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(genes.in.go.enrichment[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
      }
      
      ## Add links to GO ids
      go.result.table.with.links[["GO ID"]] <- sapply(X = go.result.table.with.links[["GO ID"]], FUN = go.link)
      
      ## Introductory text for GO enrichment table
      go.table.text <- "The table below summarizes the result of the GO term
      enrichment analysis. Each row represents a GO term significantly enriched in the selected
      gene set. The first column represents the GO term
      identifier. The second column contains a human readable description. For more details on the
      corresponding GO term, click on the identifier in the first column. The third and fourth
      column presents the p-value and q-value (adjusted p-value or FDR) capturing the level
      of significance. The fifth column displays the corresponding enrichment value E (m/n; M/N) where
      n is the number of genes with annotation from the target set, N is the number of genes with
      annotation from the gene universe, m is the number of genes from the target set annotated with the
      corresponding GO term and M is the number of genes from the gene universe annotated with
      the GO term associated with the corresponding row. The enrichment is then computed as
      E = (m/n) / (M/N). Finally, the last column, contains the genes from the target set
      annotated with the GO term represented in the corresponding row."
      
      output$textGOTable <- renderText(expr = go.table.text)
      
      ## Output table with GO enrichment result
      output$output_go_table <- renderDataTable({
        ## Error message for the user
        validate(
          need(input$selected.multiple.tfs, "Please select some transcription factor")
        )
        go.result.table.with.links #go.result.table
      },escape=FALSE,options =list(pageLength = 5)) 
      
      output$download_ui_for_go_table<- renderUI(
        tagList(downloadButton(outputId= "downloadGOData", "Download GO Enrichment"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
      )
      
      output$downloadGOData<- downloadHandler(
        filename= function() {
          paste0(paste(c("GO_enrichment",input$selected.tfs),collapse = "_"), ".tsv")
        },
        content= function(file) {
          write.table(go.result.table, 
                      file=file, 
                      sep = "\t", 
                      quote = FALSE,
                      row.names = FALSE)
        })
      
      ## Download result
      output$downloadData<- downloadHandler(
        filename= function() {
          paste("godata" , ".csv", sep="")
        },
        content= function(file) {
          write.csv(go.result.table,
                    file,
                    row.names=TRUE
          )
        })
      
      ## Link to REVIGO 
      revigo.data <- paste(revigo.data <- apply(go.result.table[,c("GO ID", "q-value")], 1, paste, collapse = " "), collapse="\n")
      
      url1 <- tags$a("here", href="#", onclick="document.revigoForm.submit();")
      url2 <- tags$form(
        name="revigoForm", action="http://revigo.irb.hr/", method="post", target="_blank",
        tags$textarea(name="inputGoList", rows="1", cols="8", class="revigoText",
                      style="visibility: hidden", revigo.data)
      )
      
      output$revigo<- renderUI(
        tagList("The enriched GO terms above may be redundant. Visualize these results in REViGO in order to remove redundancy. Click", url1, url2)
      )
      
      
      output$barplot_text <- renderText("In the following barplot each bar represents a significantly enriched 
GO term. The length of the bar corresponds to the number of genes in the
                                        target set annotated with the given GO term. The bar color captures the level
                                        of significance from blue, less significant, to red, more significant.")
      
      
      ## GO map
      output$gomap_text <- renderText("The following figure corresponds to a subgraph
                                      induced by most significant GO terms from 
                                      Biological Process subcategory. Enriched terms 
                                      are colored and the color depends on the 
                                      adjusted p-value according to the Benjamini & Hochberg 
                                      method, increasing the p-value from purple to red")
      
      output$gomap <- renderPlot(
        width     = 1040,
        height    = 1000,
        res       = 120,
        expr = {
          ## Error message for the user
          validate(
            need(input$selected.multiple.tfs, "Please select some transcription factor")
          )
          goplot(enrich.go,showCategory = 10)
        })
      
      ## Barplot
       output$bar.plot <- renderPlot(
        width     = 870,
        height    = 600,
        res       = 120,
        expr = {
          ## Error message for the user
          validate(
            need(input$selected.multiple.tfs, "Please select some transcription factor")
          )
          barplot(enrich.go,drop=TRUE,showCategory = 10)
        })
       
       ##EMAP plot
       output$emapplot_text <- renderText("The following figure consists of an enrichment map where nodes represent enriched GO terms. The
        size of a node is proportional to the number of genes annotated with the corresponding GO term in the target set.
The node colors represent the level of significance from less signficant in blue to more significant in red. Edges are drawn
between two nodes when the corresponding GO terms are semantically related.")
       
       output$emap.plot <- renderPlot(
         width     = 870,
         height    = 600,
         res       = 120,
         expr = {
           emapplot(enrich.go)
         })

       ## CNET plot
       output$cnetplot_text <- renderText("The following figure corresponds to a gene-concept network. The beige
nodes represents GO terms and the grey nodes genes. An edge is drawn from a gene to a GO term when the gene is annotated
with the corresponding gene. The size of nodes representing GO terms is proportional to the number of genes annotated
with the corresponding GO term.")

       output$cnet.plot <- renderPlot(
         width     = 870,
         height    = 600,
         res       = 120,
         expr = {
           cnetplot(enrich.go)
         })
    }
  })
  
  ##Perform KEGG pathway enrichment analysis when button is clicked
  observeEvent(input$pathway_button,{
    ## Sanity checks
    validate(
      need(input$selected.multiple.tfs, "Please select a set of transcription factors")
    )
    
    ## Determine targets of selected TFs
    # selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    # selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,input$selected.multiple.tfs] != 0)
    
    if(length(input$selected.multiple.tfs) > 1)
    {
      #gene.selection <- rowSums(selected.tfs.adj) == length(input$selected.multiple.tfs)
      
      if(input$selection.mode == "Common gene targets")
      {
        gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.agi)
      } else
      {
        gene.selection <- rowSums(selected.tfs.adj) >= 1
      }

    } else if (length(input$selected.multiple.tfs) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]    
    
    ## Set the background to perform pathway enrichment analysis depending on the user selection
    if (input$pathway_background == "allgenome")
    {
      pathway.universe <- ath.universe
    } else
    {
      pathway.universe <- network.data$name
    }
    
    ## Show element when kegg pathway enrichment starts
    shinyjs::showElement(id = 'loading.div.kegg')

    ## Compute KEGG pathway enrichment
    pathway.enrichment <- enrichKEGG(gene = selected.genes.df$name, 
                                     organism = "ath", 
                                     keyType = "kegg",
                                     universe = pathway.universe,
                                     qvalueCutoff = 0.05)
    pathway.enrichment.result <- as.data.frame(pathway.enrichment)
    
    ## Hide loading element when KEGG pathway enrichment is done
    shinyjs::hideElement(id = 'loading.div.kegg')
    
    ## Generate output table
    if(nrow(pathway.enrichment.result) > 0)
    {
      pathways.enrichment <- compute.enrichments(gene.ratios = pathway.enrichment.result$GeneRatio,
                                                 bg.ratios = pathway.enrichment.result$BgRatio)
      
      ## Separate genes with blank spaces
      kegg.enriched.genes <- pathway.enrichment.result$geneID
      for(i in 1:length(kegg.enriched.genes))
      {
        kegg.enriched.genes[i] <- paste(strsplit(kegg.enriched.genes[i],split="/")[[1]],collapse=" ")
      }
      
      ## Generate data frame with output table
      pathways.result.table <- data.frame(pathway.enrichment.result$ID, pathway.enrichment.result$Description,
                                          pathway.enrichment.result$pvalue, pathway.enrichment.result$qvalue,
                                          pathways.enrichment, 
                                          kegg.enriched.genes,
                                          stringsAsFactors = FALSE)
      
      colnames(pathways.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                           "Enrichment (Target Ratio; BG Ration)","Genes")
      
      kegg.result.table.with.links <- pathways.result.table
      
      ## Add links to genes
      for(i in 1:length(kegg.enriched.genes))
      {
        kegg.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(kegg.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
      }
      
      ## Add links to kegg pathways
      kegg.result.table.with.links[["KEGG ID"]] <- sapply(X=kegg.result.table.with.links[["KEGG ID"]],FUN = kegg.pathway.link)

      output$output_pathway_table <- renderDataTable({
        kegg.result.table.with.links
      },escape=FALSE,options =list(pageLength = 5)) 
      
      
      output$download_ui_for_kegg_table<- renderUI(
        tagList(downloadButton(outputId= "downloadKEGGData", "Download KEGG Pathway Enrichment"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
      )
      
      output$downloadKEGGData<- downloadHandler(
        filename= function() {
          paste0(paste(c("KEGG_enrichment",input$selected.tfs),collapse = "_"), ".tsv")
        },
        content= function(file) {
          write.table(pathways.result.table, 
                      file=file, 
                      sep = "\t", 
                      quote = FALSE,
                      row.names = FALSE)
        })
      
      ## Download result
      output$downloadData<- downloadHandler(
        filename= function() {
          paste("keggdata" , ".csv", sep="")
        },
        content= function(file) {
          write.csv(pathways.result.table,
                    file,
                    row.names=TRUE
          )
        })
    } else
    {
      output$no_kegg_enrichment <- renderText(expr = "No enriched KEGG pathway was detected in the selected genes.")
    }

    ## Visualization of specific enriched pathways
    genes.pathway <- rep(0, length(pathway.universe))
    names(genes.pathway) <- pathway.universe
    
    genes.pathway[selected.genes.df$name] <- 1
    
    if(nrow(pathway.enrichment.result) > 0)
    {
      pathways.for.select <- paste(pathways.result.table[["KEGG ID"]], pathways.result.table[["Description"]], sep=" - ")
      
      output$kegg_selectize <- renderUI({
        selectInput(inputId = "kegg_pathway", 
                    label="Choose Pathway for Representation",
                    multiple = FALSE,
                    selected = pathways.for.select[1],
                    choices=pathways.for.select)
      })
    ## Enriched pathway image
      output$kegg_image <- renderImage({
        pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
                 pathway.id = strsplit(input$kegg_pathway,split=" - ")[[1]][1],
                 species = "ath",
                 limit = list(gene=max(abs(genes.pathway)), cpd=1),
                 gene.idtype ="kegg")
         
        list(src = paste(c(strsplit(input$kegg_pathway,split=" - ")[[1]][1],"pathview","png"), collapse="."),
             contentType="image/png",width=1200,height=900)
      },deleteFile = T)
    }
    

    ## KEGG module enrichment analysis
    modules.enrichment <- enrichMKEGG(gene = selected.genes.df$name, 
                                      universe = pathway.universe, 
                                      organism = "ath", 
                                      keyType = "kegg",
                                      minGSSize = 4)
    
    modules.enrichment.result <- as.data.frame(modules.enrichment)
    if(nrow(modules.enrichment.result) > 0)
    {
      modules.enrichment <- compute.enrichments(gene.ratios = modules.enrichment.result$GeneRatio,
                                                bg.ratios = modules.enrichment.result$BgRatio)
      
      modules.enriched.genes <- modules.enrichment.result$geneID
      for(i in 1:length(modules.enriched.genes))
      {
        modules.enriched.genes[i] <- paste(strsplit(modules.enriched.genes[i],split="/")[[1]],collapse=" ")
      }
      
      modules.result.table <- data.frame(modules.enrichment.result$ID, modules.enrichment.result$Description,
                                         modules.enrichment.result$pvalue, modules.enrichment.result$qvalue,
                                         modules.enrichment, 
                                         modules.enriched.genes,
                                         stringsAsFactors = FALSE)
      
      colnames(modules.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                          "Enrichment (Target Ratio; BG Ration)","Genes")
      
      modules.result.table.with.links <- modules.result.table
      
      ## Add links to genes
      for(i in 1:length(modules.enriched.genes))
      {
        modules.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(modules.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
      }
      
      ## Add links to kegg pathways
      modules.result.table.with.links[["KEGG ID"]] <- sapply(X=modules.result.table.with.links[["KEGG ID"]],FUN = kegg.module.link)
      
      ## Generate output table
      output$output_module_table <- renderDataTable({
        modules.result.table.with.links
      },escape=FALSE,options =list(pageLength = 5)) 
    } else
    {
      output$text_module_kegg <- renderText(expr = "No enriched KEGG module was detected in the selected genes.")
    }
      
    
    

    
  })
  
  ##Perform TFBS enrichment analysis when button is clicked
  observeEvent(input$tfbs_button,{
    ## Sanity checks
    validate(
      need(input$selected.multiple.tfs, "Please select a set of transcription factors")
    )
    
    ## Determine targets of selected TFs
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    if(length(selected.tfs.with.zts) > 1)
    {
      
      if(input$selection.mode == "Common gene targets")
      {
        gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.agi)
      } else
      {
        gene.selection <- rowSums(selected.tfs.adj) >= 1
      }

      #gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.with.zts)
    } else if (length(selected.tfs.with.zts) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]  
    
    if(input$peak != "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, (peak.zt == input$peak & trough.zt == input$trough))
    } else if(input$peak == "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, trough.zt == input$trough)
    } else if(input$peak != "any" && input$trough == "any")
    {
      selected.genes.df <- subset(selected.genes.df, peak.zt == input$peak)
    }
    
    ## Set the background to perform TFBS enrichment analysis depending on the user selection
    if (input$tfbs_background == "allgenome")
    {
      tfbs.universe <- 33323
    } else
    {
      tfbs.universe <- 5778
    }
    
    # ## Show element when TFBS enrichment starts
    # shinyjs::showElement(id = 'loading.div.tfbs')
    
    file.precomputed <- paste0(c("precomputed_",input$up_promoter,"_",
                                 input$down_promoter, "_",
                                 input$score, "_",tfbs.universe,".tsv"),collapse="")
    
    ## Load file with precomputed results (background) and compute m and n
    precomputed.result <- read.table(file=paste0("data/precomputed_results_tfbs/",file.precomputed),header = T)
    m <- colSums(precomputed.result > 0) 
    n <- nrow(precomputed.result) - m
    
    ## Compute selection size (k) and number of ocurrences (x).
    # target.genes <- read.table(file = "peak_ZT0_trough_ZT12.txt",as.is = T)[[1]]
    target.genes <- intersect(rownames(precomputed.result),selected.genes.df$names)
    
    
    k <- length(target.genes)
    x <- colSums(precomputed.result[target.genes,] > 0)
    
    ## Compute p-values for enrichment according to a hypergeometric distribution
    p.values <- vector(mode="numeric", length=length(x))
    names(p.values) <- colnames(precomputed.result)
    
    for(i in 1:length(x))
    {
      p.values[i] <- phyper(q = x[i] - 1, m = m[i], n = n[i], k = k, lower.tail = F)
    }
    
    which(p.values < input$motif_significance)
    p.values[which(p.values < input$motif_significance)]
    
    ## Adjust p-values using Benjamini Hochberg
    q.values <- p.adjust(p = p.values,method = "BH")
    
    which(q.values < input$motif_significance)
    q.values[which(q.values < input$motif_significance)]
    
    ## Compute enrichments
    enrichments <- (x / k) / (m / nrow(precomputed.result))
    
    ## Final motifs, pvalues, qvalues and enrichments
    input <- list(motif_significance = 0.05, enrichment_threshold = 2 )
    
    sig.enrich.motifs <- names(which(q.values < input$motif_significance & enrichments > input$enrichment_threshold))
    sig.enrich.ids <- motif.ids[sig.enrich.motifs]
    final.q.values <- q.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
    final.p.values <- p.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
    final.enrichments <- enrichments[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
    
    ## Determine genes for each motif
    genes.with.motif <- vector(length = length(sig.enrich.motifs))
    for (i in 1:length(sig.enrich.motifs))
    {
      print(i)
      rows.with.motif <- which(precomputed.result[,sig.enrich.motifs[i]] != 0)
      all.genes.with.motif <- rownames(precomputed.result)[rows.with.motif]
      genes.with.motif[i] <- paste(... = intersect(all.genes.with.motif,target.genes), collapse = ",")
    }
    
    ## Motifs logos
    motifs.images <- paste0("motifs_images/",sig.enrich.motifs)
    
    for (i in 1:length(motifs.images))
    {
      motifs.images[i] <- paste0("<img src='",motifs.images[i],".png', align = 'center', width = 100>")
    }
    
    
    ## Store data
    tfbs.result.table <- data.frame(sig.enrich.motifs, sig.enrich.ids, motifs.images, final.p.values, final.q.values, final.enrichments, genes.with.motif, stringsAsFactors = FALSE) 
    colnames(tfbs.result.table) <- c("DNA motifs", "Motif ID", "DNA logo", "P-values", "Q-values", "Enrichments", "Genes")
    
    ## Add links to jaspar motifs
    tfbs.result.table[["Motif ID"]] <- sapply(X=sig.enrich.ids,FUN = tfbs.link)
    
    ## Add links to genes
    for (i in 1:length(genes.with.motif))
    {
      tfbs.result.table$Genes[i] <- paste(sapply(X = strsplit(genes.with.motif[i], split = ",")[[1]],FUN = gene.link.function), collapse = ", ")
    }
    
    tfbs.result.table <- tfbs.result.table[order(final.q.values),]
    ## Output table with TFBS enrichment result
    output$output_tfbs_table <- renderDataTable({
      tfbs.result.table
    },escape = FALSE,options =list(pageLength = 10))
    
    output$download_ui_tfbs_table<- renderUI(
      tagList(downloadButton(outputId= "downloadTFBSData", "Download TFBS Enrichment"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
    )
    
    output$downloadTFBSData<- downloadHandler(
      filename= function() {
        paste0(paste(c("TFBS_enrichment",selected.tfs.with.zts, input$peak, input$trough),collapse = "_"), ".tsv")
      },
      content= function(file) {
        write.table(tfbs.result.table, 
                    file=file, 
                    sep = "\t", 
                    quote = FALSE,
                    row.names = FALSE)
      })
    
    
    
    
  })

}

# Run the application 
shinyApp(ui = ui, server = server)

