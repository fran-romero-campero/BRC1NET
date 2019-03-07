## Load libraries
library(shiny)
library(ggplot2)
library(org.At.tair.db)
library(SuperExactTest)

##Auxiliary functions
intersect2sets <- function(set1, set2, alias, gene.descriptions){
  intersection.data <- list()
  sets <- list(set1, set2)
  results <- supertest(x = sets, n = 5778)
  results.table <- summary(results)
  p.value <- tail(results.table$P.value, n=1) #Get the last p-value
  enrichment <- (results.table$Table)[["FE"]][nrow(results.table$Table)]
  intersection.genes <- (results.table$Table)[["Elements"]][nrow(results.table$Table)]
  intersection.genes <- strsplit(intersection.genes, split = ", ")[[1]]
  
  intersection.genes.agi <- intersection.genes
  intersection.genes.primary.symbol <- alias[intersection.genes]
  names(intersection.genes.primary.symbol) <- NULL
  gene.table <- matrix(nrow=length(intersection.genes), ncol=3)
  colnames(gene.table) <- c("AGI", "Primary Symbol", "Description")
  gene.table[,1] <- intersection.genes.agi
  gene.table[,2] <- intersection.genes.primary.symbol
  #  gene.table[,3] <- description
  
  
  
  
  intersection.genes.description <- gene.descriptions[intersection.genes]
  names(intersection.genes.description) <- NULL
  
  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment
  intersection.data[[3]] <- data.frame(intersection.genes,intersection.genes.primary.symbol,intersection.genes.description,stringsAsFactors = F)
  
  names(intersection.data) <- c("p-value", "enrichment", "gene.table")
  return(intersection.data)
  
}

## Load network
network.data <- read.table(file="data/brc1_transcriptional_network.tsv",header = TRUE,as.is=TRUE,sep="\t",comment.char = "")
gene.names <- network.data$names

## Extract gene descriptions
gene.description <- network.data$T.annotation
names(gene.description) <- network.data$names

#head(network.data)
colnames(network.data)[56] <- "ATHB-53"
colnames(network.data)[54] <- "DOF5-4"
colnames(network.data)[51] <- "NAC6-ORE1"
colnames(network.data)[46] <- "ATHB-40"
colnames(network.data)[34] <- "ATHB-6"
colnames(network.data)[33] <- "ATHB-21"

## Transforming coordinates for a better visualization
x.coord <- as.numeric(network.data$y.pos)
y.coord <- as.numeric(network.data$x.pos)

## Extract gene ids
genes <- sort(network.data$name)

## Create tair links for genes
gene.links <- vector(mode="character",length=length(genes))
for(i in 1:length(genes))
{
  tair.link <- paste0("https://www.arabidopsis.org/servlets/TairObject?type=locus&name=",genes[i])
  gene.links[i] <- paste(c("<a href=\"",
                       tair.link,
                       "\" target=\"_blank\">",
                       genes[i], "</a>"),
                     collapse="")
}


names(gene.links) <- genes

## Load all and circadian genes
my.key <- keys(org.At.tair.db, keytype="ENTREZID")
my.col <- c("SYMBOL", "TAIR")
columns(org.At.tair.db)
alias2symbol.table <- select(org.At.tair.db, keys=my.key, columns=my.col, keytype="ENTREZID")
alias2symbol.table <- subset(alias2symbol.table, genes %in% TAIR)
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias[is.na(alias)] <- "" 
genes.selectize <- paste(names(alias), alias, sep=" - ")

## Extract TFs
tfs.data <- read.table(file = "data/network_tfs.txt",header=T, as.is=T)
tf.ids <- tfs.data$AGI
names(tf.ids) <- tfs.data$name

#For loop to classify all genes in clusters
brc1.clusters <- rep( list(list()), 9 ) 

for (i in 1:(nrow(network.data)-1))
{
  print(i)
  if (strsplit(x = network.data$T.cluster[i], split = "|")[[1]][1] == 1)
  {
    brc1.clusters[[1]] <- append(brc1.clusters[[1]], gene.names[i])
  }
  
  if (strsplit(x = network.data$T.cluster[i], split = "|")[[1]][6] == 2)
  {
    brc1.clusters[[2]] <- append(brc1.clusters[[2]], gene.names[i])
  }
  
  if (strsplit(x = network.data$T.cluster[i], split = "|")[[1]][11] == 3)
  {
    brc1.clusters[[3]] <- append(brc1.clusters[[3]], gene.names[i])
  }
  
  if (strsplit(x = network.data$T.cluster[i], split = "|")[[1]][16] == 4)
  {
    brc1.clusters[[4]] <- append(brc1.clusters[[4]], gene.names[i])
  }
  
  if (strsplit(x = network.data$T.cluster[i], split = "|")[[1]][21] == 5)
  {
    brc1.clusters[[5]] <- append(brc1.clusters[[5]], gene.names[i])
  }
  
  if (strsplit(x = network.data$T.cluster[i], split = "|")[[1]][26] == 6)
  {
    brc1.clusters[[6]] <- append(brc1.clusters[[6]], gene.names[i])
  }
  
  if (strsplit(x = network.data$T.cluster[i], split = "|")[[1]][31] == 7)
  {
    brc1.clusters[[7]] <- append(brc1.clusters[[7]], gene.names[i])
  }
  
  if (strsplit(x = network.data$T.cluster[i], split = "|")[[1]][36] == 8)
  {
    brc1.clusters[[8]] <- append(brc1.clusters[[8]], gene.names[i])
  }
  
  if (strsplit(x = network.data$T.cluster[i], split = "|")[[1]][41] == 9)
  {
    brc1.clusters[[9]] <- append(brc1.clusters[[9]], gene.names[i])
  }
}


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
  titlePanel(tags$b("A Transcriptional Network Downstream of BRC1")),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        
        tags$h3(tags$b("Gene Target Selection:")),
        
        checkboxGroupInput(inputId = "selected.tfs",
                           label = "Select Transcription Factors:",
                           choices = sort(tfs.data$name,decreasing=F),
                           inline = TRUE,width = "100%"),
        
        checkboxInput(inputId = "edges",label = "Visualize edges",value = FALSE),
        
        actionButton(inputId = "button_tfs",label = "Select Genes"),
        
        selectizeInput(inputId = "topological_parameter",
                       label = "Choose topological parameter",
                       choices = c("In degree",
                                   "Transitivity",
                                   "Closeness",
                                   "Betweeness",
                                   "Eccentricity")),
        
        conditionalPanel(condition = "input.topological_parameter == 'In degree'",
                         sliderInput(inputId = "indegree_range",
                                     label="Choose an in degree range:",
                                     min=min(network.data$indegree),
                                     max=max(network.data$indegree),
                                     value=c(min(network.data$indegree),max(network.data$indegree))),
                         actionButton(inputId = "button_indegree",label = "Select Genes")
                         ),

        conditionalPanel(condition = "input.topological_parameter == 'Transitivity'",
                         sliderInput(inputId = "transitivity_range",
                                     label="Choose a transitivity range:",
                                     min=min(network.data$transitivity),
                                     max=max(network.data$transitivity),
                                     value=c(min(network.data$transitivity),max(network.data$transitivity))),
                         actionButton(inputId = "button_transitivity",label = "Select Genes")
        ),

        conditionalPanel(condition = "input.topological_parameter == 'Closeness'",
                         sliderInput(inputId = "closeness_range",
                                     label="Choose a closeness range:",
                                     min=min(network.data$closeness),
                                     max=max(network.data$closeness),
                                     value=c(min(network.data$closeness),max(network.data$closeness))),
                         actionButton(inputId = "button_closeness",label = "Select Genes")
        ),
        
        conditionalPanel(condition = "input.topological_parameter == 'Betweeness'",
                         sliderInput(inputId = "betweeness_range",
                                     label="Choose a betweeness range:",
                                     min=min(network.data$betweenness),
                                     max=max(network.data$betweenness),
                                     value=c(min(network.data$betweeness),max(network.data$betweeness))),
                         actionButton(inputId = "button_betweeness",label = "Select Genes")
        ),
        

        conditionalPanel(condition = "input.topological_parameter == 'Eccentricity'",
                         sliderInput(inputId = "eccentricity_range",
                                     label="Choose an eccentricity range:",
                                     min=min(network.data$eccentricity),
                                     max=max(network.data$eccentricity),
                                     value=c(min(network.data$eccentricity),max(network.data$eccentricity))),
                         actionButton(inputId = "button_eccentricity",
                                      label = "Select Genes")
        ),
        
        tags$h3(tags$b("Intersections")),
        selectInput(inputId = "cluster", label="Cluster",
                    choices = 1:9, selected = NULL,
                    multiple = FALSE, selectize = TRUE), 
        selectInput(inputId = "top_parameter", label = "Topological parameter", 
                    choices = c("Degree","Betweeness", "Closeness", "Eccentricity","Transitivity"), selected = NULL,
                    multiple = FALSE, selectize = TRUE),
        selectInput(inputId = "threshold", label = "Parameter threshold",
                    choices = c(0.75,0.90,0.95), selected = NULL, multiple = FALSE, selectize = TRUE),
        
        
        actionButton(inputId = "top_intersect", label = "Test"),
        
                
                
#         sliderInput(inputId = "degree_range", label = h3("Degree Range"), min = 0, 
#                     max = 11, value = c(2, 4)),
#         actionButton(inputId = "button_degree",label = "Select Genes"),


        width = 3 
      ),
      
      # Show a plot of the network and table with selected genes info
      mainPanel(
         plotOutput("networkPlot",hover="plot_hover"),#click="plot_click"),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         tags$br(),
         dataTableOutput(outputId = "output_table"),
         
         width = 9
      )
   )
)

# Define server logic required to represent the network
server <- function(input, output) {

  ## Default network representation 
  default.representation <- ggplot(network.data, aes(x.pos,y.pos)) + 
    theme(panel.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.y = element_blank()) + 
    geom_point(color=network.data$color,size=1)
 
  ## Initial/default visualization of BRC1 Network
  output$networkPlot <- renderPlot({
    print("default")
    default.representation
  },height = 700)
  
  
#  observeEvent(input$plot_click, {
    output$networkPlot <- renderPlot({
      print("click")
      default.representation +
        annotate("text",
                 x = input$plot_hover$x + 1,
                 y = input$plot_hover$y + 1, 
                 label="lala",size=10)
    },height = 700)
  
  #})  

  ## Visualization of selected genes according to their degree
  observeEvent(input$button_tfs, {
    print("aquí llego 4")
    

    if(length(input$selected.tfs) == 1)
    {
      gene.selection <- network.data[,input$selected.tfs] == 1
      sum(gene.selection)
    } else if(length(input$selected.tfs) > 1)
    {
      gene.selection <- rowSums(network.data[,input$selected.tfs]) == length(input$selected.tfs)
      
    } else if(length(input$selected.tfs) == 0)
    {
      gene.selection <- NULL
    }
    
    selected.genes.df <- network.data[gene.selection,]
    selected.nodes.colors <- selected.genes.df$color
    
    selected.tfs.df <- subset(network.data, names %in% tf.ids[input$selected.tfs])
    
    if(input$edges)
    {
      network.representation <- ggplot(network.data, aes(x.pos,y.pos)) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank())
      
      network.representation <- ggplot(network.data, aes(x.pos,y.pos)) + 
          theme(panel.background = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks.y = element_blank()) + 
          geom_point(color=network.data$color,size=1) +
          geom_point(data = selected.tfs.df, size=8, fill=selected.tfs.df$color,colour="black",pch=21) +
          geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21) 
      
      for(i in 1:length(input$selected.tfs))
      {
        tf.xpos <- subset(network.data, names == tf.ids[input$selected.tfs[i]])[["x.pos"]]
        tf.ypos <- subset(network.data, names == tf.ids[input$selected.tfs[i]])[["y.pos"]]
        network.representation <- network.representation +
            annotate("segment",
                     x=rep(tf.xpos,nrow(selected.genes.df)),
                     y=rep(tf.ypos,nrow(selected.genes.df)),
                     xend=selected.genes.df$x.pos,
                     yend=selected.genes.df$y.pos, 
                     color="grey", arrow=arrow(type="closed",length=unit(0.1, "cm")))
      }
      
      output$networkPlot <- renderPlot({
        network.representation
      },height = 700)
    } else
    {
      output$networkPlot <- renderPlot({
        ggplot(network.data, aes(x.pos,y.pos)) + 
          theme(panel.background = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks.y = element_blank()) + 
          geom_point(color=network.data$color,size=1) +
          geom_point(data = selected.tfs.df, size=8, fill=selected.tfs.df$color,colour="black",pch=21) +
          geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21) 
      },height = 700)
    }

    
    
    selected.genes.df$names <- gene.links[selected.genes.df$names]
    
    output$output_table <- renderDataTable({
    #  selected.genes.df[,c("names","S.name","S.annotation","S.TF/Other","T.cluster")]#as.data.frame(genes.annotation.data.with.links)
      selected.genes.df[,c("names","S.name","S.annotation","T.mapman","T.TF.Other","T.cluster")]
    },escape=FALSE)  
  })
  

  
  
  ## Visualization of selected genes according to their degree
  observeEvent(input$button_indegree, {
    print("aquí llego 5")

    degree.values <- input$indegree_range
    selected.genes.df <- subset(network.data, indegree >= degree.values[1] & indegree <= degree.values[2])
    selected.nodes.colors <- selected.genes.df$color

    output$networkPlot <- renderPlot({
      ggplot(network.data, aes(x.pos,y.pos)) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank()) +
        geom_point(color=network.data$color,size=1) +
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    },height = 700)


    selected.genes.df$names <- gene.links[selected.genes.df$names]

    output$output_table <- renderDataTable({
      #  selected.genes.df[,c("names","S.name","S.annotation","S.TF/Other","T.cluster")]#as.data.frame(genes.annotation.data.with.links)
      selected.genes.df[,c("names","S.name","S.annotation","T.mapman","T.TF.Other","T.cluster")]
    },escape=FALSE)
  })
  
  ## Visualization of selected genes according to their transitivity
  observeEvent(input$button_transitivity, {
    print("aquí llego 6")
    
    trans.values <- input$transitivity_range
    selected.genes.df <- subset(network.data, transitivity >= trans.values[1] & transitivity <= trans.values[2])
    selected.nodes.colors <- selected.genes.df$color
    
    output$networkPlot <- renderPlot({
      ggplot(network.data, aes(x.pos,y.pos)) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank()) +
        geom_point(color=network.data$color,size=1) +
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    },height = 700)
    
    
    selected.genes.df$names <- gene.links[selected.genes.df$names]
    
    output$output_table <- renderDataTable({
      #  selected.genes.df[,c("names","S.name","S.annotation","S.TF/Other","T.cluster")]#as.data.frame(genes.annotation.data.with.links)
      selected.genes.df[,c("names","S.name","S.annotation","T.mapman","T.TF.Other","T.cluster")]
    },escape=FALSE)
  })
  
  
  ## Visualization of selected genes according to their closeness
  observeEvent(input$button_closeness, {
    print("aquí llego 7")
    
    close.values <- input$closeness_range
    selected.genes.df <- subset(network.data, closeness >= close.values[1] & closeness <= close.values[2])
    selected.nodes.colors <- selected.genes.df$color
    
    output$networkPlot <- renderPlot({
      ggplot(network.data, aes(x.pos,y.pos)) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank()) +
        geom_point(color=network.data$color,size=1) +
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    },height = 700)
    
    
    selected.genes.df$names <- gene.links[selected.genes.df$names]
    
    output$output_table <- renderDataTable({
      #  selected.genes.df[,c("names","S.name","S.annotation","S.TF/Other","T.cluster")]#as.data.frame(genes.annotation.data.with.links)
      selected.genes.df[,c("names","S.name","S.annotation","T.mapman","T.TF.Other","T.cluster")]
    },escape=FALSE)
  })
  
  ## Visualization of selected genes according to their betweeness
  observeEvent(input$button_betweeness, {
    print("aquí llego 8")
    
    betw.values <- input$betweeness_range
    selected.genes.df <- subset(network.data, betweenness >= betw.values[1] & betweenness <= betw.values[2])
    selected.nodes.colors <- selected.genes.df$color
    output$networkPlot <- renderPlot({
      ggplot(network.data, aes(x.pos,y.pos)) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank()) +
        geom_point(color=network.data$color,size=1) +
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    },height = 700)
    
    
    selected.genes.df$names <- gene.links[selected.genes.df$names]
    
    output$output_table <- renderDataTable({
      #  selected.genes.df[,c("names","S.name","S.annotation","S.TF/Other","T.cluster")]#as.data.frame(genes.annotation.data.with.links)
      selected.genes.df[,c("names","S.name","S.annotation","T.mapman","T.TF.Other","T.cluster")]
    },escape=FALSE)
  })

  
  ## Visualization of selected genes according to their betweeness
  observeEvent(input$button_eccentricity, {
    print("aquí llego 8")
    
    eccen.values <- input$eccentricity_range
    selected.genes.df <- subset(network.data, eccentricity >= eccen.values[1] & eccentricity <= eccen.values[2])
    selected.nodes.colors <- selected.genes.df$color
    output$networkPlot <- renderPlot({
      ggplot(network.data, aes(x.pos,y.pos)) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank()) +
        geom_point(color=network.data$color,size=1) +
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    },height = 700)
    
    
    selected.genes.df$names <- gene.links[selected.genes.df$names]
    
    output$output_table <- renderDataTable({
      #  selected.genes.df[,c("names","S.name","S.annotation","S.TF/Other","T.cluster")]#as.data.frame(genes.annotation.data.with.links)
      selected.genes.df[,c("names","S.name","S.annotation","T.mapman","T.TF.Other","T.cluster")]
    },escape=FALSE)
  })
  
  ##Visualization and intersection between topological parameters and gene clusters
  observeEvent(input$top_intersect, {
    print("Test top intersection")
    gene.names <- network.data$names
    
    if (input$top_parameter == "Degree")
    {
      brc1net.degree <- network.data$indegree + network.data$outdegree
      degree.threshold <- quantile(brc1net.degree, prob=as.numeric(input$threshold))
      top.genes <- gene.names[brc1net.degree > degree.threshold]
    } else if (input$top_parameter == "Transitivity")
    {
      network.data$transitivity[is.na(network.data$transitivity)] <- 0
      brc1net.trans <- network.data$transitivity
      trans.threshold <- quantile(brc1net.trans, prob=as.numeric(input$threshold))
      top.genes <- gene.names[brc1net.trans > trans.threshold]
    } else if (input$top_parameter == "Closeness")
    {
      brc1net.closeness <- network.data$closeness
      closeness.threshold <- quantile(brc1net.closeness, prob=as.numeric(input$threshold))
      top.genes <- gene.names[brc1net.closeness > closeness.threshold]
    } else if (input$top_parameter == "Betweeness")
    {
      brc1net.bet <- network.data$betweeness
      bet.threshold <- quantile(brc1net.bet, prob=as.numeric(input$threshold))
      top.genes <- gene.names[brc1net.bet > bet.threshold]
    } else if (input$top_parameter == "Eccentricity")
    {
      brc1net.eccen <- network.data$eccentricity
      eccen.threshold <- quantile(brc1net.eccen, prob=as.numeric(input$threshold))
      top.genes <- gene.names[brc1net.eccen > eccen.threshold]
    } 
    
    if (input$cluster == 1)
    {
      cluster <- brc1.clusters[[1]]
    } else if (input$cluster == 2)
    {
      cluster <- brc1.clusters[[2]]
    } else if (input$cluster == 3)
    {
      cluster <- brc1.clusters[[3]]  
    } else if (input$cluster == 4)
    {
      cluster <- brc1.clusters[[4]]  
    } else if (input$cluster == 5)
    {
      cluster <- brc1.clusters[[5]]  
    } else if (input$cluster == 6)
    {
      cluster <- brc1.clusters[[6]]  
    } else if (input$cluster == 7)
    {
      cluster <- brc1.clusters[[7]]  
    } else if (input$cluster == 8)
    {
      cluster <- brc1.clusters[[8]]  
    } else if (input$cluster == 9)
    {
      cluster <- brc1.clusters[[9]]  
    }
    
    result <- intersect2sets(set1 = top.genes, set2 = cluster, alias = alias, gene.descriptions = gene.description)
    p.value <- result[1][[1]]
    enrichment <- result[2][[1]]
    intersect.genes <- result[3][[1]]$intersection.genes
    
    
    
    selected.genes.df <- subset(network.data, names %in% intersect.genes)
    # selected.nodes.colors <- selected.colors[selected.genes.df$peak.zt]
    
    # selected.genes.df <- network.data[gene.selection,]
    selected.nodes.colors <- selected.genes.df$color
    
    
    output$networkPlot <- renderPlot({
      ggplot(network.data, aes(x.pos,y.pos)) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank()) +
        geom_point(color=network.data$color,size=1) +
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    },height = 700)
    
    
    ## Visualization of text with p value and enrichment
    if (length(top.genes) == 0)
    {
      text.intersection.result <- "<b>There is no genes with this restriction<b>"
    } else 
    {
      if(p.value < 0.01)
      {
        text.intersection.result <- paste0("<b>The intersection between the genes with high ", input$topological_parameter,
                                           " and genes that show a ", input$peak_top, " peak and ", input$trough_top, " trough ",
                                           " is significant with a p-value of ", p.value,
                                           " and an enrichment of ", round(x = enrichment,digits = 2),
                                           "<b>") 
        
      } else
      {
        text.intersection.result <- paste0("<b>The intersection between the genes with high ", input$topological_parameter,
                                           " and genes that show a ", input$peak_top, " peak and ", input$trough_top, " trough ",
                                           " is NOT significant with a p-value of ", round(x=p.value, digits = 2),
                                           " and an enrichment of ", round(x = enrichment,digits = 2),
                                           "<b> <br> <br>") 
      }
      
    }
    
    
    output$outputText <- renderText(expr = text.intersection.result, quoted = FALSE)
    
    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

