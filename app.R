## Load libraries
library(shiny)
library(ggplot2)
library(org.At.tair.db)


## Load network
network.data <- read.table(file="data/brc1_transcriptional_network.tsv",header = TRUE,as.is=TRUE,sep="\t",comment.char = "")
head(network.data)


## Tranforming coordinates for a better visualization
x.coord <- as.numeric(network.data$y.pos)
y.coord <- as.numeric(network.data$x.pos)

## Extract gene ids
genes <- sort(network.data$name)

## Create tair links for genes

i <- 1
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
        
        # conditionalPanel(condition = "input.topological_parameter == 'Transitivity'",
        #                  sliderInput(inputId = "transitivity_range",
        #                              label="Choose a transitivity range:",
        #                              min=min(network.data$transitivity),
        #                              max=max(network.data$transitivity),
        #                              value=c(min(network.data$transitivity),max(network.data$transitivity))),
        #                  actionButton(inputId = "button_transitivity",label = "Select Genes")
        # ),

#         sliderInput(inputId = "degree_range", label = h3("Degree Range"), min = 0, 
#                     max = 11, value = c(2, 4)),
#         actionButton(inputId = "button_degree",label = "Select Genes"),


        width = 3 
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("networkPlot"),
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

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ## Initial/default visualization of ATTRACTOR
  output$networkPlot <- renderPlot({
    ggplot(network.data, aes(x.pos,y.pos)) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank()) + 
      geom_point(color=network.data$color,size=1)
  },height = 700)
  

  ## Visualization of selected genes according to their degree
  observeEvent(input$button_degree, {
    print("aquí llego 3")
    
    node.degree <- network.data$indegree
    degree.values <- input$degree_range
    selected.genes.df <- subset(network.data, indegree >= degree.values[1] & indegree <= degree.values[2])
    selected.nodes.colors <- selected.genes.df$color
    
    print(selected.genes.df)
    print(selected.nodes.colors)
    
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
  })

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
      
    }
    selected.genes.df <- network.data[gene.selection,]
    selected.nodes.colors <- selected.genes.df$color
    
    selected.tfs.df <- subset(network.data, names %in% tf.ids[input$selected.tfs])

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
  
  
  
  
    
  
}

# Run the application 
shinyApp(ui = ui, server = server)

