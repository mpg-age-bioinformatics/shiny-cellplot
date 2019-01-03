.libPaths("/srv/shiny-server/cellplot/libs")
gitversion <- function(){ 
  git<-read.csv("/srv/shiny-server/.git/refs/heads/master", header=FALSE)
  git<-git$V1
  git<-toString(git[1])
  git<-substr(git, 1, 7)
  return(git)
}
library(shiny)
library(CellPlot)
futile.logger::flog.threshold(futile.logger::ERROR, name = "cellplotLogger")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  # reformat input data
  plot.data <- reactive({
    inFile <- input$file1
    filetype <- input$filetype
    req(inFile)
    req(input$file2)
    
    filetype_map <- c("xlsx" = 'xlsx',  'tsv' = '\t', 'csv' = ',', 'txt'=" ")
    if(filetype == 'auto'){
      file_extension =  unlist(strsplit(inFile$datapath, '[.]'))[length(unlist(strsplit(inFile$datapath, '[.]')))]
      if(file_extension %in% names(filetype_map)){
        filetype <- filetype_map[file_extension]
        names(filetype) <- NULL
        
      } else {
        print(paste("wrong file format", file_extension))
        return(NULL)
      }
    }
    
    if(filetype == 'xlsx'){
      D <- read.xlsx(inFile$datapath, sheetIndex = 1, header = input$header)
    } else {
      D <- read.csv(inFile$datapath, header = input$header, sep = filetype)
    }
    
    Categories<-unique(D$Category)
    
    inFileExp <- input$file2
    filetype2 <- input$filetype2

    filetype_map <- c("xlsx" = 'xlsx',  'tsv' = '\t', 'csv' = ',', 'txt'=" ")
    if(filetype2 == 'auto'){
      file_extension =  unlist(strsplit(inFileExp$datapath, '[.]'))[length(unlist(strsplit(inFileExp$datapath, '[.]')))]
      if(file_extension %in% names(filetype_map)){
        filetype2 <- filetype_map[file_extension]
        names(filetype2) <- NULL
        
      } else {
        print(paste("wrong file format", file_extension))
        return(NULL)
      }
    }
    
    if(filetype2 == 'xlsx'){
      DD <- read.xlsx(inFileExp$datapath, sheetIndex = 1, header = input$header2)
    } else {
      DD <- read.csv(inFileExp$datapath, header = input$header2, sep = filetype2)
    }
    
    DD[DD == ''] <- NA

    vars <- names(DD)
    if (is.null(input$categories)){
      updateSelectInput(session, "categories","Select Categories", choices = Categories, selected = NULL)
    }
    if ( input$genessel == ""){
      updateSelectInput(session, "genessel","Select Genes Name/ID Column", choices = vars, selected = NULL )
    }
    if (input$logfcsel == ""){
      updateSelectInput(session, "logfcsel","Select Log2(FC) Column", choices = vars, selected = NULL)
    }
    if (input$padjsel == ""){
      updateSelectInput(session, "padjsel","Select P Adj. Column", choices = vars, selected = NULL)
    }
    
    req(!is.null(input$categories))
    req( input$genessel != "") 
    req( input$logfcsel != "") 
    req( input$padjsel != "")
    
    siggenes<-DD[input$genessel]
    DD['GenesSignificant']<-siggenes
    #dput(DD,"/srv/shiny-server/cellplot/test.ref.R")
    
    D<-D[D["Category"] == input$categories, ]
    D <- D[order(D$PValue),]
    D <- D[1:input$nterms,]
    
    D["GenesSignificant"]<-D$Genes
    genes<-D$GenesSignificant
    class(genes)
    genes<-lapply( genes, function(x) strsplit(toString(x), ", ")[[1]])
    D$GenesSignificant<-genes
    
    D$log2FoldChange<-lapply( genes, function(x) DD[DD$GenesSignificant %in% x, input$logfcsel ])
    D$padj<- lapply( genes, function(x) DD[DD$GenesSignificant %in% x, input$padjsel ])
    
    D$LogEnrich<-D$Fold.Enrichment
    
    return(D)
  })
  
  plot.cellplot<-reactive({
    
    x<-plot.data()
    cell.plot(x = setNames(x$LogEnrich, x$Term), 
              cells = x$log2FoldChange, 
              main ="GO enrichment", 
              x.mar = c(.5, 0.1),
              key.n = 7, 
              y.mar = c(.1, 0), 
              cex = 1.6, 
              cell.outer = 3, 
              bar.scale = .7, 
              space = .2)
  })
  
  plot.symplot<-reactive({
    
    x<-plot.data()
    sym.plot(x = setNames(x$LogEnrich, x$Term), 
             cells = x$log2FoldChange, 
             x.annotated = x$Count, 
             main = "GO enrichment",
             x.mar = c(.7, 0.1),
             y.mar = c(0.1,0.1),
             key.n = 7, 
             cex = 1.6, 
             axis.cex = .8, 
             group.cex = .7) 
  })
  
  plot.arclot<-reactive({
    
    x<-plot.data()
    x$up <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i>0] })
    x$dwn <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i<0] })
    arc.plot(x = setNames(x$LogEnrich, x$Term), 
             up.list = x$up, 
             down.list = x$dwn, 
             x.mar = c(0.7, 0.1), # c(0.7, 0.1)
             y.mar = c(0.3, 0.1)) # c(0.1, 0.1)
    
  })
  
  #plot.histogram<-reactive({
  #  x<-plot.data()
  #  y <- lapply(x, function (x) {
  #    x$Upregulated <- sapply(x$log2FoldChange, function (z) sum(z>0))
  #    x$Downregulated <- sapply(x$log2FoldChange, function (z) sum(z<0))
  #    x
  #  })
  #  yterms <- unique(unlist(lapply(y, function(x){
  #    x <- subset(x, pvalCutOff <= 0.05)
  #    x <- x[order(x$LogEnrich),]
  #    head(x, 9)$GO.ID
  #  })))
  #  
  #  par(mar = c(0,.5,2.5,8))
  #  go.histogram(y, go.alpha.term = "pvalCutOff", gene.alpha.term = "padj", 
  #               min.genes = 5, max.genes = 1e10, go.selection = yterms, show.ttest = T,
  #               main = "GO enrichment", 
  #               axis.cex = 1, lab.cex = 1.5, main.cex = 1.5)
  #})

  output$cellplot <- renderPlot({
    plot.cellplot()
  })
  
  output$downloadcellPlot <- downloadHandler(
    # specify file name
    filename = function(){
      paste0(input$outfile,".cellplot.",gitversion(),'.pdf')
    },
    content = function(filename){
      pdf(filename,height = 12.75, width = 13.50)
      x<-plot.data()
      cell.plot(x = setNames(x$LogEnrich, x$Term), 
                cells = x$log2FoldChange, 
                main ="GO enrichment", 
                x.mar = c(.5, 0.1),
                key.n = 7, 
                y.mar = c(.1, 0.1), 
                cex = 1.6, 
                cell.outer = 3, 
                bar.scale = .7, 
                space = .2)
      dev.off()
    }
    
  )
  
  output$symplot <- renderPlot({
    plot.symplot()
  }) 
  
  output$downloadsymPlot <- downloadHandler(
    # specify file name
    filename = function(){
      paste0(input$outfile,".symplot.",gitversion(),'.pdf')
    },
    content = function(filename){
      pdf(filename,height = 12.75, width = 13.50)
      x<-plot.data()
      sym.plot(x = setNames(x$LogEnrich, x$Term), 
               cells = x$log2FoldChange, 
               x.annotated = x$Count, 
               main = "GO enrichment",
               x.mar = c(.7, 0.1),
               y.mar = c(0.1,0.1),
               key.n = 7, 
               cex = 1.6, 
               axis.cex = .8, 
               group.cex = .7) 
      dev.off()
    }  
    )
  
  
  output$arcplot <- renderPlot({
    plot.arclot()
  })
  
  output$downloadarcPlot <- downloadHandler(
    # specify file name
    filename = function(){
      paste0(input$outfile,".arcplot.",gitversion(),'.pdf')
    },
    content = function(filename){
      pdf(filename,height = 12.75, width = 13.50)
      x<-plot.data()
      x$up <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i>0] })
      x$dwn <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i<0] })
      arc.plot(x = setNames(x$LogEnrich, x$Term), 
               up.list = x$up, 
               down.list = x$dwn, 
               x.mar = c(0.7, 0.1), # c(0.7, 0.1)
               y.mar = c(0.3, 0.1)) # c(0.1, 0.1)
      dev.off()
    }  
  )
  
  #output$histogram <- renderPlot({
  #  plot.histogram()
  #})
  
  output$appversion <- renderText ({ 
    paste0('App version: <b>',gitversion(),'</b>')
  }
  )
})