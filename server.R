.libPaths("/srv/shiny-server/cellplot/libs")
gitversion <- function(){ 
  git<-read.csv("/srv/shiny-server/.git/modules/cellplot/refs/heads/master", header=FALSE)
  git<-git$V1
  git<-toString(git[1])
  git<-substr(git, 1, 7)
  return(git)
}
library(shiny)
library(CellPlot)
futile.logger::flog.threshold(futile.logger::ERROR, name = "cellplotLogger")

cell.plot_ = function(
  x, cells, x.col=NULL, cell.col=c("deepskyblue2","white","coral"),
  inf.shading = 30/cell.lwd,  space=0.1, x.mar=c(0.2, 0), y.mar = c(0.1, 0), x.bound=NULL, lab.cex = 1, xdes.cex=1, xlab.cex=1, xlab.ticks=5,
  sym=FALSE, cell.lwd=2, cell.outer=2, cell.sort=T, cell.limit=30, cell.bounds=NULL, elem.bounds=NULL, xlab="GO Term Enrichment",
  key=T, key.lab="Differential Expression", key.n=11, spacers=NULL, bar.scale=1, gridlines=T, cellborder=NA, both.sides = 'yes', ... )
{
  # parameter checks
  if(!is.null(x.bound)){ if(!(is.numeric(x.bound) && (x.bound > 0)) ) {
    stop("x.bound must be a positive numeric value")
  }}
  
  if(!is.null(elem.bounds)) {
    ec = sapply(cells,length)
    excl = which( (ec < elem.bounds[1]) | (ec > elem.bounds[2]) )
    if (length(excl) == length(x) ) { stop("No elements in the specified range!") }
    if (length(excl) > 0) {
      x = x[-excl]
      if (!is.null(cells)) { cells = cells[-excl] } 
      else { x.up=x.up[-excl]; x.down=x.down[-excl] }
      if (!is.null(x.col)) { x.col = x.col[-excl] }
    }
  }
  
  # yscale = par("pin")[1]/par("pin")[2] * diff(par("usr")[1:2])/diff(par("usr")[3:4])
  yscale = (diff(par("usr")[3:4])/par("pin")[2])
  
  par(xpd=NA)
  cell.col.inf = cell.col[c(1,3)]
  #if (is.null(xlab.ticks)) { xlab.ticks = round( max(x) / 10, digits = 1) }
  #ticksize = xlab.ticks
  ybound = c(1,0) + c(-1,1)*y.mar
  
  # scale  
  ysteps = ybound[1] - cumsum( rep(0.3 * yscale * bar.scale, length(x)+1+ifelse(is.null(spacers), 0, length(spacers)) ) )
  ybound[2] = min(ysteps, ybound[2])
  if ( ybound[2] < par("usr")[3] ) {
    warning("Plotting area too small! Decrease bar.scale.fixed or increase vertical space.")
  }
  
  xbound = c(0,1) + c(1,-1)*x.mar
  if (is.null(spacers)) {
    ysteps = seq( ybound[1], ybound[2], length.out=( length(x)+1 ) )
  } else {
    spacers = spacers + 1:length(spacers) + 1
    ysteps = seq( ybound[1], ybound[2], length.out=( length(x)+1+length(spacers) ) )
    ysteps = ysteps[-spacers]
  }
  
  ygap = abs(ysteps[1]-ysteps[2])
  yspace = space * ygap
  
  celldata = unlist(cells)
  cellinf <- is.infinite(celldata)
  cellmis <- is.na(celldata)
  cellbound = range(celldata[!cellinf & !cellmis])
  if (!is.null(cell.bounds) & is.numeric(cell.bounds) ) { cellbound = cell.bounds }
  if (sym) { cellbound = rep( max(abs(cellbound)),2 ) * c(-1,1) }
  # cellcolmap = seq( cellbound[1], cellbound[2], length.out=101 )  
  
  # inlcude here, if color should be centered around 0 or only up/down is being visualized
  
  if (both.sides == 'no'){
    cellcolmap = c(seq(cellbound[1], cellbound[2], length.out = 50))
    names(cellcolmap) = c(colorRampPalette(c(cell.col[1], cell.col[2]))(50))
    
    
  } else {
    cellcolmap = c( seq( cellbound[1], 0, length.out=50 ), 0, seq( 0, cellbound[2], length.out=50 ) )
    names(cellcolmap) = c( colorRampPalette( c(cell.col[1], cell.col[2]) )(50),
                           cell.col[2],
                           colorRampPalette( c(cell.col[2], cell.col[3]) )(50) )
    
  }

  if (any(cellinf)) {
    if (key.n < 5) stop("key.n must be >4 when infinite values are present")
  } else {
    if (key.n < 3) stop("key.n must be >2")
  }
  
  
  labvec = rep("",length(x))
  if (!is.null(names(x))) { labvec=names(x) }
  colvec = x.col
  if( is.null(colvec) ) { colvec = rep("black",length(x)) }
  
  # do the actual plotting
  plot.new()
  
  axis.at <- seq(xbound[1], xbound[2], length.out = xlab.ticks)
  axis.lab <- round(seq(min(0,min(x)), max(x), length.out = xlab.ticks),1)
  if (!is.null(x.bound) && is.numeric(x.bound) && (x.bound > 0)) { axis.lab <- round(seq(min(0,min(x)), x.bound, length.out = xlab.ticks),1) }
  
  # GRID
  if (gridlines) {
    segments(axis.at, ybound[2]-yspace, axis.at, ybound[1]+0.015, col="grey", lwd = cell.outer )
    #segments(left.axis.at[length(left.axis.at)],ybound[2]-yspace,left.axis.at[1],ybound[2]-yspace, col="grey", lwd = bar.lwd)
    segments(axis.at[length(axis.at)],ybound[2]-yspace,axis.at[1],ybound[2]-yspace, col="grey", lwd = cell.outer)
  }
  axis(3, pos=ybound[1]+0.015, at = axis.at, labels = axis.lab, cex.axis=xlab.cex, padj=NA, lwd=cell.outer )
  
  for (i in 1:length(x)) {
    bar.n <- length(cells[[i]])
    bar.nreal <- sum(!is.na(cells[[i]]))
    xsteps = seq(xbound[1], (x[i]/max(x))*(xbound[2]-xbound[1])+xbound[1], length.out=(bar.n+1))
    if (!is.null(x.bound) && is.numeric(x.bound) && (x.bound > 0)) { 
      xsteps = seq(xbound[1], (x[i]/x.bound)*(xbound[2]-xbound[1])+xbound[1], length.out=(bar.n+1)) }
    
    if (is.null(cells)) { xsteps = c(xbound[1], (x[i]/max(x))*(xbound[2]-xbound[1])+xbound[1]) }
    # row labels
    text( xbound[1], ysteps[i+1]+ygap*0.5, labels=labvec[i], pos=2, cex=lab.cex, col = colvec[i] )
    if(cell.sort) { cells[[i]] = sort(cells[[i]]) }
    bar.order <- order(cells[[i]])
    # number of cells labels
    text( xsteps[length(xsteps)], ysteps[i+1]+ygap*0.5, pos=4, cex=lab.cex, labels=bar.nreal)
    # map colors to cell values
    bar.val <- cells[[i]]
    if (bar.nreal < 1) bar.val <- 0
    bar.inf <- is.infinite(bar.val)
    bar.inf.pos <- bar.inf & bar.val > 0
    bar.inf.neg <- bar.inf & bar.val < 0
    bar.i <- sapply(bar.val, function(bi) which.min(abs(cellcolmap - bi)))
    bar.col <- names(cellcolmap)[bar.i]
    #bar.shade <- ifelse(bar.inf, 10L, NA_integer_)
    
    #bar.cell.lwd <- ifelse( bar.n < cell.limit, cell.lwd, NA)
    bar.cell.lwd = cell.lwd
    bar.shade <- ifelse(bar.inf, inf.shading, 0)
    bar.col.inf = bar.col
    if (bar.n < cell.limit) { bar.col.inf = "black" } else { bar.col.inf[bar.inf] = "black" }
    bar.col[bar.inf.neg] <- cell.col.inf[1]
    bar.col[bar.inf.pos] <- cell.col.inf[2]
    
    # omit cell borders if the bar has more than cell.limit cells
    # bar.border <- if (bar.n < cell.limit) rep("black",bar.n) else bar.col
    # make little boxes
    rect(xsteps[-(bar.n+1)], ysteps[i+1]+ygap-yspace, xsteps[-1], ysteps[i+1]+yspace,
         col = bar.col, lwd = bar.cell.lwd, border = NA) # works only whith lwd >0
    rect(xsteps[-(bar.n+1)], ysteps[i+1]+ygap-yspace, xsteps[-1], ysteps[i+1]+yspace,
         col = bar.col.inf, lwd = bar.cell.lwd, density = bar.shade, border = cellborder) # works only whith lwd >0
    # and another box around the whole bar
    rect( xbound[1], ysteps[i+1]+ygap-yspace, xsteps[length(xsteps)], ysteps[i+1]+yspace, lwd=cell.outer )
  }
  
  
  
  title( ... )
  # XAXIS DESIGNATION
  text( (xbound[1]+xbound[2])/2, ybound[1]+0.015+strheight("0",cex = xlab.cex)*2, labels=xlab, pos=3, cex=xdes.cex )
  segments( xbound[1], ybound[1]+0.015, xbound[1], ybound[2]-yspace, lwd=cell.outer)
  
  # color legend
  if (key) {
    lc = c( 0.8,ybound[2]-ygap-yspace,0,ybound[2]-yspace*3 )
    absmax = max(abs(cellbound))
    lc.min <- min(cellbound)
    lc.max <- max(cellbound)
    cellcolmap.center <- cellcolmap[51]
    cellcolmap <- cellcolmap[lc.min <= cellcolmap & cellcolmap <= lc.max]
    lc.xsteps = seq( xbound[1], xbound[2], length.out=key.n+1 )
    lc.xgap = lc.xsteps[1] - lc.xsteps[2]
    
    lc.density <- rep(0, key.n)
    if (any(cellinf)) {
      lc.range <- c(-Inf, seq( lc.min, lc.max, length.out=key.n-2), Inf)
      lc.col <- c(cell.col.inf[1], names(cellcolmap)[seq(1,length(cellcolmap),length.out=key.n-2)], cell.col.inf[2])
      lc.density[c(1,length(lc.density))] = inf.shading
      if (key.n %% 2 > 0 && lc.min < 0 && lc.max > 0) {
        lc.range = c( -Inf, seq( lc.min, 0, length.out=(key.n-1)/2 ), seq( 0, lc.max, length.out=(key.n-1)/2 )[-1], Inf )
        print(lc.range)
      }
    } else {
      lc.range = seq( lc.min, lc.max, length.out=key.n )
      lc.col <- names(cellcolmap)[seq(1,length(cellcolmap),length.out=key.n)]
      if (key.n %% 2 > 0 && lc.min < 0 && lc.max > 0) {
        lc.range = c( seq( lc.min, 0, length.out=(key.n-1)/2+1 ), seq( 0, lc.max, length.out=(key.n-1)/2+1 )[-1] )
        lc.col <- names(cellcolmap)[seq(1,length(cellcolmap),length.out=key.n)]
      }
    }
    
    key.num = as.character(round(lc.range,1))
    if (!is.null(cell.bounds)) {
      or = range(celldata[!cellinf & !cellmis])
      if (cell.bounds[1] > min(or)) { key.num[1+(any(cellinf))] = paste0("<",key.num[1+(any(cellinf))])}
      if (cell.bounds[2] < max(or)) { key.num[length(key.num)-(any(cellinf))] = paste0(key.num[length(key.num)-(any(cellinf))],">")}
    }
    
    rect( lc.xsteps[-(key.n+1)]-lc.xgap*.1, lc[2], lc.xsteps[-1]+lc.xgap*.1, lc[4], col=lc.col, lwd=cell.outer )
    rect( lc.xsteps[-(key.n+1)]-lc.xgap*.1, lc[2], lc.xsteps[-1]+lc.xgap*.1, lc[4], col="black", lwd=cell.lwd, density = lc.density, border=NA )
    text( (lc.xsteps[-(key.n+1)]+lc.xsteps[-1])/2, lc[2], pos=1, labels=key.num, cex=xlab.cex, font=2 )
    text( (xbound[1]+xbound[2])/2, lc[2]-strheight("0",cex=xlab.cex)*1.5 , labels=key.lab, pos=1, cex=xdes.cex )
  }
  par(xpd=T)
}





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
    #if (input$padjsel == ""){
    #  updateSelectInput(session, "padjsel","Select P Adj. Column", choices = vars, selected = NULL)
    #}
    
    # need to select a category otherwise no plot
    validate(
      need(!is.null(input$categories), paste("please select at least one category"))
    )
    
    req(!is.null(input$categories))
    req( input$genessel != "") 
    req( input$logfcsel != "") 
    #req( input$padjsel != "")
    
    siggenes<-DD[input$genessel]
    DD['GenesSignificant']<-siggenes
    #dput(DD,"/srv/shiny-server/cellplot/test.ref.R")
    
    #cat(file=stderr(), input$categories, "\n")
    D<-D[D$Category %in% input$categories, ]
    #D <- D[order(D$PValue),]
    nterms=input$nterms
    if (nterms >= nrow(D)) nterms = nrow(D)
    D <- D[1:nterms,]
    
    D["GenesSignificant"]<-D$Genes
    genes<-D$GenesSignificant
    class(genes)
    genes<-lapply( genes, function(x) strsplit(toString(x), ", ")[[1]])
    D$GenesSignificant<-genes
    
    D$log2FoldChange<-lapply( genes, function(x) DD[DD$GenesSignificant %in% x, input$logfcsel ])
    #D$padj<- lapply( genes, function(x) DD[DD$GenesSignificant %in% x, input$padjsel ])
    
    D$LogEnrich<-D$Fold.Enrichment
    
    return(D)
  })
  
  # color
  cell.col <- reactive({
    cell.col = gsub('[ ]', '', input$cell.col)
    cell.col = unlist(strsplit(cell.col, ','))  
  })

  plot.cellplot<-reactive({

    x<-plot.data()

    x_mar = 1 - input$x_mar
    y_mar = 1 - input$y_mar
    cell.plot_(x = setNames(x$LogEnrich, x$Term), 
              cells = x$log2FoldChange, 
              main = input$text.main, 
              x.mar = c(x_mar, 0),
              key.n = 7, 
              y.mar = c(y_mar, 0), 
              lab.cex = input$label_size,
              xlab.cex = input$label_size,
              xdes.cex = input$label_size,
              key.lab = input$text.key,
              xlab = input$text.xaxis,
              both.sides = input$both.sides,
              cell.col = cell.col(),
              cex = 1.6, 
              cell.outer = 2, 
              bar.scale = .7, 
              space = .2,
              cell.limit = input$cell.limit,
              cell.lwd = input$cellbordersize,
              cellborder=input$cellborder)
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
      x_mar = 1 - input$x_mar
      y_mar = 1 - input$y_mar
      cell.plot_(x = setNames(x$LogEnrich, x$Term), 
                 cells = x$log2FoldChange, 
                 main = input$text.main, 
                 x.mar = c(x_mar, 0),
                 key.n = 7, 
                 y.mar = c(y_mar, 0), 
                 lab.cex = input$label_size,
                 xlab.cex = input$label_size,
                 xdes.cex = input$label_size,
                 key.lab = input$text.key,
                 xlab = input$text.xaxis,
                 both.sides = input$both.sides,
                 cell.col = cell.col(),
                 cex = 1.6, 
                 cell.outer = 2, 
                 bar.scale = .7, 
                 space = .2,
                 cell.limit = input$cell.limit,
                 cell.lwd = input$cellbordersize,
                 cellborder=input$cellborder)
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
