library(conflicted)
library(shiny)
library(ggplot2)
library(plotly)
library(shinycssloaders)
library(plyr)
library(heatmaply)
library(threejs)

mychoices = c("0","1","2","3","4","5","6","7")

heatmapInfo <- conditionalPanel(condition = "input.inTabset==1", h4("HEATMAPINFO"))
vulcanoInfo <- conditionalPanel(condition = "input.inTabset==3", h4("VULCANOINFO"))

heatNetPane2 <- sidebarPanel(
  conditionalPanel(condition = "input.inTabset==2 || input.inTabset == 1",
  fluidRow(
    column(6,checkboxInput("NAWM", "NAWM")),
    column(6,selectInput("NAWM_REG", NULL, choices = c("UP", "DOWN", "BOTH"), selected = "BOTH"))
  ),
  fluidRow(
    column(6, checkboxInput("AL", "active", value = T)),
    column(6, selectInput("AL_REG", NULL, choices = c("UP", "DOWN", "BOTH"), selected = "UP"))
  ),
  
  fluidRow(
    column(6, checkboxInput("IL", "inactive", value = T)),
    column(6, selectInput("IL_REG", NULL, choices = c("UP", "DOWN", "BOTH"), selected = "DOWN"))
  ),
  
  fluidRow(
    column(6,checkboxInput("CA", "chronic active", value = T)),
    column(6,selectInput("CA_REG", NULL, choices = c("UP", "DOWN", "BOTH"), selected = "DOWN"))
  ),
  
  fluidRow(
    column(6,checkboxInput("RL", "remyelinating")),
    column(6,selectInput("RL_REG", NULL, choices = c("UP", "DOWN", "BOTH"), selected = "BOTH"))
  ),
  tags$hr(style="border-color: black; margin-top: 0.1em; margin-bottom: 0.5em"),
  fluidRow(
    column(6, textInput("FDR", "FDR < ", value = "0.05")),
    column(6, radioButtons("FDR_ALL", label = "",
                           choices = c("at least one", "all")))
  ),
  fluidRow(
    column(6, textInput("logFC", "abs(logFC) >", value = "0")),
    column(6, radioButtons("LOGFC_ALL", label = "",
                           choices = c("at least one", "all")))
  ),
  tags$hr(style="border-color: black; margin-top: 0.1em; margin-bottom: 0.5em"),
  fluidRow(
    column(12, checkboxInput("adjust_colors", "Adjust Color Range")),
    column(12, checkboxInput("row_dendogram", "Dendogram")),
    column(12, textInput("n_genes", "N most significant genes", value = "inf"))
  ),
  fluidRow(
    conditionalPanel(condition = "input.inTabset==2",
    column(12, tags$hr(style="border-color: black; margin-top: 0.1em; margin-bottom: 0.5em")),
    
    column(12, radioButtons(inputId = "enrichment", label = "Network Enrichment",
                          choices = c("Connected Component", "KeyPathwayMiner"))),
    
    column(6, conditionalPanel(condition="input.enrichment == 'Connected Component'",
                               selectInput("C", "Component", choices = mychoices, selected = "1"))
    ),
    
    column(6, conditionalPanel(condition = "input.enrichment == 'KeyPathwayMiner'",
                               selectInput("K", "K", choices = mychoices, selected = "0"))),
    
    column(12,actionButton ("applyToNetwork", HTML("<b>","Apply Changes", "</b>"), 
                                         style="color: white; background-color: #b70101; border-color: black; width: 100%"))
    
  )),
  tags$hr(style="border-color: black; margin-top: 0.1em; margin-bottom: 0.5em"),
  
  fluidRow(
    textAreaInput("geneInput", "Gene Symbol")
  ),
  
  fluidRow(
    column(12,
           fileInput("geneList", "Choose a CSV File", multiple = F, accept = "text/csv")
    )
  ),
  tags$hr(style="margin-top: 0.1em; margin-bottom: 0.5em"),
  fluidRow(
    column(12,
           actionButton("start", HTML("<b>","Start", "</b>"), 
                        style="color: white; background-color: #b70101; border-color: black; width: 100%")
    )
  ),
  tags$hr(style="border-color: black; margin-bottom: 0.5em"),
  fluidRow(
    column(6,
           downloadButton("export_genes", "Genes")
    ),
    
    column(6,
           downloadButton("export_network", "Network")
    )
  )),
  conditionalPanel(condition = "input.inTabset==3",
                   radioButtons("volcanoSlection", label = "Lesion Type",
                                choices = c("NAWM", "active", "inactive", "chronic active", "remyelinating"),
                                selected = "NAWM")),
  width = 2
)

ui <- navbarPage("Multiple Sclerosis Gene Atlas", id = "navID", theme = "style.css",
                 tabPanel("Home", id="tabHome",
                          includeMarkdown("home.md")), 
                 tabPanel("Get Started", id = "tabID",
                          # Application title
                          #titlePanel("Multiple Sklerosis Genome Atlas"),
                          sidebarLayout(
                            #pane1,
                            heatNetPane2,
                            
                            mainPanel(id = "mainPanel", heigth = "100%",
                              tabsetPanel(type = "pills",
                                tabPanel("HeatMap", value = "1", plotlyOutput("heatMap", height = "750px") %>% withSpinner()),
                                tabPanel("Network", value = "2", scatterplotThreeOutput("network", height = "750px") %>% withSpinner()),
                                tabPanel("Volcano", value = "3", plotlyOutput("volcano", height = "750px") %>% withSpinner()),
                                id = "inTabset",
                                selected = "1"
                              ),
                              width = 10
                            )
                          )
                 ),
                 
                 tabPanel("Guide",
                          includeMarkdown("guide.md")),
                 
                 tabPanel("About",
                          includeMarkdown("about.md")
                          )
)
