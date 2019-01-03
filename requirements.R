.libPaths("/srv/shiny-server/cellplot/libs")
library(shiny)
# apt-get install libssl-dev

if(!require(devtools)){
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}

if(!require(CellPlot)){
  library(devtools)
  devtools::install_github("dieterich-lab/CellPlot", build_vignettes = TRUE)
  library(CellPlot)
}

if(!require(futile.logger)){
  install.packages("futile.logger", dependencies = TRUE)
  library(futile.logger)
}


quit(save="no")


#################################################
#												
#												
#                   DEV AREA					
#												
#												
#################################################


library(CellPlot)

data("leukemiasGO")

head(leukemiasGO$CLL)



test.data<-leukemiasGO$CLL
write.table(test.data,"/srv/shiny-server/cellplot/test.data.txt")
str(test.data)
# GenesAnnotated
# GenesSignificant
# padj
# log2FoldChange

#test.data$GenesSignificant # list
#typeof(sig.genes)
#typeof(sig.genes[702][1])


#######################################
sig.genes<-test.data$GenesSignificant
sig.genes<-lapply( sig.genes, function(x) toString(x))
sig.genes<- unlist(sig.genes, use.names=FALSE)

padj.genes<-test.data$padj
padj.genes<-lapply( padj.genes, function(x) toString(x))
padj.genes<- unlist(padj.genes, use.names=FALSE)

log2FoldChange.genes<-test.data$log2FoldChange
log2FoldChange.genes<-lapply( log2FoldChange.genes, function(x) toString(x))
log2FoldChange.genes<- unlist(log2FoldChange.genes, use.names=FALSE)

df<-data.frame(GenesSignificant=sig.genes,padj=padj.genes,log2FoldChange=log2FoldChange.genes)
write.table(df, "/srv/shiny-server/cellplot/genes.padj.logfc.txt", sep="\t")
#######################################

ref_table<-read.csv("/srv/shiny-server/cellplot/ref.values.txt",sep="\t")
head(ref_table)

#######################################

df<-data.frame(GO.ID=test.data$GO.ID,Term=test.data$Term,Annotated=test.data$Annotated,Significant=test.data$Significant,
               Expected=test.data$Expected,pvalCutOff=test.data$pvalCutOff,LogEnrich=test.data$LogEnrich,GenesSignificant=sig.genes)
write.table(tt,"/srv/shiny-server/cellplot/test.input.DAVID.txt",sep="\t")
df<-read.table("/srv/shiny-server/cellplot/test.input.DAVID.txt")
#######################################

#genes<-tt$GenesSignificant
##genes<-genes[as.list(strsplit(toString(genes), ", ")[[1]])]
#class(genes)
##genes<-list(genes)
#genes
##genes<-lapply( genes, function(x) as.list(strsplit(toString(x), ", ")[[1]]) )
#genes<-lapply( genes, function(x) strsplit(toString(x), ", ")[[1]])
#genes
#typeof(genes[702][1])

#######################################

library(CellPlot)

#######################################

df<-read.table("/srv/shiny-server/cellplot/test.input.DAVID.txt")
df<-read.table("/srv/shiny-server/cellplot/DAVIDws.CLL.cellplot.FunctionalAnnotation.398ecca.tsv", sep="\t", header=TRUE)
df$GenesSignificant<-df$Genes
genes<-df$GenesSignificant
class(genes)
genes<-lapply( genes, function(x) strsplit(toString(x), ", ")[[1]])
df$GenesSignificant<-genes 
str(df)

#######################################

ref_table<-read.csv("/srv/shiny-server/cellplot/ref.values.tsv",sep="\t")
head(ref_table)
df$log2FoldChange<- lapply( genes, function(x) ref_table[ref_table$GenesSignificant %in% x, 'log2FoldChange'])
df$padj<- lapply( genes, function(x) ref_table[ref_table$GenesSignificant %in% x, 'padj'])
df$LogEnrich<-df$Fold.Enrichment

#######################################
dput(x,file="/srv/shiny-server/cellplot/test.object.R")
testx<-dget("/srv/shiny-server/cellplot/test.object.R")
View(testx)

testr<-dget("/srv/shiny-server/cellplot/test.ref.R")
View(testr)


#x <- subset(df, pvalCutOff <= 0.05 & Significant > 20)

x <- df[order(df$PValue),]
x <- x[1:20,]

cell.plot(x = setNames(x$LogEnrich, x$Term), 
          cells = x$log2FoldChange, 
          main ="GO enrichment", 
          x.mar = c(.4, 0), 
          key.n = 7, 
          y.mar = c(.1, 0), 
          cex = 1.6, 
          cell.outer = 3, 
          bar.scale = .7, 
          space = .2)

sym.plot(x = setNames(x$LogEnrich, x$Term), 
         cells = x$log2FoldChange, 
         x.annotated = x$Count, 
         main = "GO enrichment",
         x.mar = c(.47, 0), 
         key.n = 7, 
         cex = 1.6, 
         axis.cex = .8, 
         group.cex = .7) 

x$up <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i>0] })
x$dwn <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i<0] })
arc.plot(x = setNames(x$LogEnrich, x$Term), 
         up.list = x$up, 
         down.list = x$dwn, 
         x.mar = c(.9, .5))



y <- lapply(leukemiasGO, function (x) {
  x$Upregulated <- sapply(x$log2FoldChange, function (z) sum(z>0))
  x$Downregulated <- sapply(x$log2FoldChange, function (z) sum(z<0))
  x
})
yterms <- unique(unlist(lapply(y, function(x){
  x <- subset(x, pvalCutOff <= 0.05)
  x <- x[order(x$LogEnrich),]
  head(x, 9)$GO.ID
})))

par(mar = c(0,.5,2.5,8))
go.histogram(y, go.alpha.term = "pvalCutOff", gene.alpha.term = "padj", 
             min.genes = 5, max.genes = 1e10, go.selection = yterms, show.ttest = T,
             main = "GO enrichment", 
             axis.cex = 1, lab.cex = 1.5, main.cex = 1.5)


tmp<-read.table("/srv/shiny-server/cellplot/DAVIDws.FunctionalAnnotation.398ecca.csv", sep="\t",header=TRUE)
Categories<-tmp$Category
Categories<-unique(Categories)
typeof(Categories)
class(Categories)
typeof(Categories[1])
typeof(unlist(Categories))







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
    helpText(a(href = "https://raw.githubusercontent.com/mpg-age-bioinformatics/shiny/master/histogram/chol.txt", "Example input")),
    selectInput("categories","Select Categories", choices = NULL ),#, multiple = TRUE),
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
    helpText(a(href = "https://raw.githubusercontent.com/mpg-age-bioinformatics/shiny/master/histogram/chol.txt", "Example input")),
    selectInput("genes","Select Genes Name/ID Column", choices = NULL),
    selectInput("log2fc","Select Log2(FC) Column", choices = NULL),
    selectInput(session, "padj","Select P Adj. Column", choices = NULL),
    hr(),
    sliderInput('nterms',  "Number of terms to plot", min = 1, max = 20, value = 10, step = 1),
    hr(),
    textInput("outfile", "Output file name", value="cellplots"),
    submitButton('generate plots')
  ),
  mainPanel(
    plotOutput("cellplot", height = "500px", width = "500px"),
    br(), br(),
    p("This App uses the", code('cellplot'), " package. For more information read the respective documentation in ",
      a("github", href = "http://htmlpreview.github.io/?https://github.com/dieterich-lab/CellPlot/blob/master/vignettes/CellPlotManual.html"),
      "."
    ),
    p("Please keep the version tag on all downloaded files."),
    htmlOutput('appversion')
  )
)






