.libPaths("/srv/shiny-server/cellplot/libs")
library(shiny)
# Define UI for application that draws a histogram
shinyUI( fluidPage(
  tags$head(
    tags$style(HTML("hr {border-top: 2px solid #000000;}"))
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "DAVID's Functional Annotation file",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",
                           ".tsv",
                           ".xlsx")
      ),
      radioButtons("filetype", "Please select DAVID's Functional Annotation file type", choices = c('auto' = 'auto', 
                                                                                                    "excel" = 'xlsx',  
                                                                                                    'tab-separated' = '\t', 
                                                                                                    'comma-seperated' = ',', 
                                                                                                    'semicolon-separated' = ';'), inline = TRUE),
      checkboxInput("header", "Header", TRUE),
      helpText(a(href = "https://raw.githubusercontent.com/mpg-age-bioinformatics/shiny/master/cellplot/DAVIDws.CLL.cellplot.FunctionalAnnotation.398ecca.tsv", "Example input")),
      hr(),
      fileInput("file2", "Log2FC and P Adj. reference file",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",
                           ".tsv",
                           ".xlsx")
      ),
      radioButtons("filetype2", "Please select Log2FC and P Adj. reference file type", choices = c('auto' = 'auto', 
                                                                                                   "excel" = 'xlsx',  
                                                                                                   'tab-separated' = '\t', 
                                                                                                   'comma-seperated' = ',', 
                                                                                                   'semicolon-separated' = ';'), inline = TRUE),
      checkboxInput("header2", "Header", TRUE),
      helpText(a(href = "https://raw.githubusercontent.com/mpg-age-bioinformatics/shiny/master/cellplot/ref.values.tsv", "Example input")),
      hr(),
      selectInput("categories","Select Categories", choices = c(NULL), multiple = TRUE, selected = NULL),
      selectInput("genessel","Select Genes Name/ID Column",  choices = c(NULL), selected = NULL),
      selectInput("logfcsel","Select Log2(FC) Column", choices = c(NULL), selected = NULL),
      selectInput("padjsel","Select P Adj. Column", choices = c(NULL), selected = NULL),
      sliderInput('nterms',  "Number of terms to plot", min = 1, max = 40, value = 40, step = 1),
      hr(),
      textInput("outfile", "Output file name", value="cellplots"),
      submitButton('submit')
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "cellplot",
          plotOutput("cellplot", height = "850px", width = "900px"),
          downloadButton('downloadcellPlot', 'Download Plot')
          ),
        tabPanel(
          "symplot",
          plotOutput("symplot", height = "1000px", width = "1000px"),
          downloadButton('downloadsymPlot', 'Download Plot')
        ),
        tabPanel(
          "arcplot",
          plotOutput("arcplot", height = "850px", width = "900px"),
          downloadButton('downloadarcPlot', 'Download Plot')
        )),
        #tabPanel(
        #  "histogram",
        #  plotOutput("histogram", height = "1000px", width = "1000px")
        #)
        #),
      br(), br(),
      p("This App uses the", code('cellplot'), " package. For more information read the respective documentation in ",
        a("github", href = "http://htmlpreview.github.io/?https://github.com/dieterich-lab/CellPlot/blob/master/vignettes/CellPlotManual.html"),
        "."
      ),
      p("Please keep the version tag on all downloaded files."),
      htmlOutput('appversion')
    )
  )
))

