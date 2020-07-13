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
      helpText(a(href = "https://datashare.mpcdf.mpg.de/s/9PLrf8rBOucc3ht/download", "Example input")),
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
      helpText(a(href = "https://datashare.mpcdf.mpg.de/s/mBp904u65knP1Fr/download", "Example input")),
      hr(),
      selectInput("xvalues","Select values for x-axis", choices = c(NULL), selected = NULL),
      selectInput("categories","Select Categories", choices = c(NULL), multiple = TRUE, selected = NULL),
      selectInput("genessel","Select Genes Name/ID Column",  choices = c(NULL), selected = NULL),
      selectInput("logfcsel","Select Log2(FC) Column", choices = c(NULL), selected = NULL),
      submitButton('submit'),
      hr(),
      #selectInput("padjsel","Select P Adj. Column", choices = c(NULL), selected = NULL),
      textInput('text.main', 'Plot title', value = "GO enrichment"),
      textInput('text.xaxis', 'X-axis title', value = "Differential Expression"),
      textInput('text.key', 'Legend title', value = "GO Term Enrichment"),
      
      sliderInput('nterms',  "Number of terms to plot", min = 1, max = 40, value = 40, step = 1),
      submitButton('submit'),
      hr(),
      textInput('cell.col', 'Input colors for cells', value = "deepskyblue2, white, coral" ),
      p("Pick colors from ", a(href = "http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf", "this table", target="_blank"), ", or input RGB values"),
      radioButtons("both.sides", "Colors centered around 0", choices = c("yes","no"), selected = "yes", inline = TRUE),
      
      sliderInput('cellbordersize',  "Border width for each cell", min = 0.5, max = 4, value = 1, step = 0.2),
      sliderInput('cell.limit', "Select cutoff for how many individual cells to plot before switching to gradient", min = 1, max = 50, step = 1, value = 10),
      checkboxInput("cellborder", "Cell border", TRUE),
      submitButton('submit'),
      hr(),
      
      sliderInput('x_mar', 'Select width of plot area', min = 0, max = 1, value = 0.35, step = 0.01),
      sliderInput('y_mar', 'Select height of plot area', min = 0, max = 1, value = 0.99, step = 0.01),
      
      sliderInput('label_size', 'Select term label size', min = 0.1, 2, value = 1, step = 0.01),
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

