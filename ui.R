## Author:  Francisco J. Romero-Campero
## Contact: Francisco J. Romero-Campero - fran@us.es 

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
                     "Multiple transcription factor analysis" = "multiple_gene",
                     "Individual gene analysis" = "individual_gene",
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
                       tags$div(align="justify", "This online tool allows researchers to explore the transcriptional
          network downstream of", tags$i("BRC1"),". This network was constructed by combining genome-wide transcriptional 
          profiling of active/dormant buds and seedlings after BRC1 induction, with 
          genome-wide BRC1 binding sites determined by Chromatin Immunoprecipitation 
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
                                BRC1-dependent transcription factors (TFs) over their common targets. The node representing BRC1 is located
                                at the center. Genes are classified into", tags$b("nine different co-expressed clusters"), "constituted by overexpressed (UP) or
                                underexpressed (DOWN) genes after BRC1 induction. Red and purple nodes are BRC1 direct targets, blue and green nodes
                                represent indirect BRC1 targets whereas purple and green nodes denote genes encoding TFs. Follow the steps below to
                                explore this network:", 
                                
                                tags$ol(
                                  tags$li("Select your TFs of interest using the", tags$b("Select Transcription Factors"),
                                          "checkbox menu on the left below."),
                                  tags$li("You can also select a gene cluster using the dropdown menu", tags$b("Select a gene cluster:")),
                                  tags$li("Check the box", tags$b("Visualize Edges"), "when you want to depict arrows from TFs to their target genes."),
                                  tags$li("Choose between two", tags$b("Modes of selction"), ", either selecting all the target genes or only those common to the selected TFs."),
                                  tags$li("Click on the", tags$b("SELECT GENES"), "to visualize your selected TFs target genes located
                                           in the specified cluster. Explore the different tabs to 
                                           download a table with the selected genes, perform a signficance analysis of the overlap between the selected 
                                           TFs targets and the specified gene cluster as well as functional enrichment analysis.")
                                )
                       )
      ),
      
      
      conditionalPanel(condition = "input.navigation_bar == 'github'",
                       tags$div(align = "justify", tags$b("BRC1NET,"), "is entirely developed using 
        the R package", tags$b( tags$a(href="https://shiny.rstudio.com/", "shiny.")), "The 
        source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                tags$b("GitHub."), "If you experience any problem using BRC1NET please create an", tags$b(tags$a(href="https://github.com/fran-romero-campero/BRC1TRANSNET/issues","issue")), "in GitHub and we will address it."),
                       tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/BRC1TRANSNET", "BRC1NET at GitHub"))))
      ),
      
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
                                                    "Cluster UP_C1" = "1", "Cluster UP_C2" = "2", "Cluster UP_C3" = "3",
                                                    "Cluster UP_C4" = "4", "Cluster UP_C5" = "5", "Cluster UP_C6" = "6",
                                                    "Cluster DOWN_C1" = "7", "Cluster DOWN_C2" = "8", "Cluster DOWN_C3" = "9"),
                                        selected = "any", multiple = FALSE, selectize = TRUE),
                            checkboxInput(inputId =  "edges",label = "Visualize Edges",value = FALSE),
                            
                            radioButtons(inputId = "selection.mode",label = "Mode of selection:",
                                         choices = c("Common gene targets","All gene targets"),selected = "Common gene targets",inline = F),
                            
                            actionButton(inputId = "go_multiple",label = "Select Genes")
                     ),
                     
                     column(width = 9,
                            tabsetPanel(type = "tabs",
                                        tabPanel(title = "Network Visualization",
                                                 plotOutput("networkPlot", click="plot_click")
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
                                                                               the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in BRC1NET:"),
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
                                                                               the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in BRC1NET:"),
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
