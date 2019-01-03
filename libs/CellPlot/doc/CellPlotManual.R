## ------------------------------------------------------------------------
library(CellPlot)
data("leukemiasGO")

## ------------------------------------------------------------------------
x <- subset(leukemiasGO$CLL, pvalCutOff <= 0.05 & Significant > 20)
x <- x[order(-x$LogEnrich),]

## ---- fig.width=7, fig.height=8------------------------------------------
cell.plot(x = setNames(x$LogEnrich, x$Term), 
          cells = x$log2FoldChange, 
          main ="GO enrichment (NoT vs CLL)", 
          x.mar = c(.4, 0), 
          key.n = 7, 
          y.mar = c(.1, 0), 
          cex = 1.6, 
          cell.outer = 3, 
          bar.scale = .7, 
          space = .2)

## ---- fig.width=7, fig.height=6------------------------------------------
sym.plot(x = setNames(x$LogEnrich, x$Term), 
         cells = x$log2FoldChange, 
         x.annotated = x$Annotated, 
         main = "GO enrichment (NoT vs CLL)",
         x.mar = c(.47, 0), 
         key.n = 7, 
         cex = 1.6, 
         axis.cex = .8, 
         group.cex = .7) 

## ---- fig.width=7, fig.height=6------------------------------------------
x$up <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i>0] })
x$dwn <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i<0] })
arc.plot(x = setNames(x$LogEnrich, x$Term), 
         up.list = x$up, 
         down.list = x$dwn, 
         x.mar = c(.9, .5))

## ------------------------------------------------------------------------
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

## ---- fig.width=7, fig.height=6------------------------------------------
par(mar = c(0,.5,2.5,8))
go.histogram(y, go.alpha.term = "pvalCutOff", gene.alpha.term = "padj", 
             min.genes = 5, max.genes = 1e10, go.selection = yterms, show.ttest = T,
             main = "GO enrichment\nin leukemia differential gene expression\ncompared to control samples", 
             axis.cex = 1, lab.cex = 1.5, main.cex = 1.5)

## ---- eval=FALSE---------------------------------------------------------
#  # CRAN
#  library(parallel)
#  library(dplyr)
#  library(stringr)
#  # Bioconductor
#  library(BiocParallel)
#  library(DESeq2)
#  library(topGO)
#  library(annotate)
#  library(leukemiasEset)
#  
#  data(leukemiasEset)
#  
#  # M <- select(hu6800.db, featureNames(leukemiasEset), c("ENSEMBL","GO"), keytype = "ENSEMBL")
#  # M <- subset(M, ONTOLOGY == "BP", c("ENSEMBL","GO"))
#  # M <- M[!duplicated(M),]
#  # M <- dlply(M, "ENSEMBL", function (x) unique(x$GO))
#  
#  # M <- as.character(unique(unlist(select(org.Hs.eg.db,featureNames(leukemiasEset),"ENSEMBL",keytype = "ENSEMBL"))))
#  
#  A <- as(leukemiasEset,"data.frame")
#  A <- subset(A, LeukemiaType %in% c("NoL","ALL", "AML", "CLL"))
#  As <- subset(A, select = "LeukemiaType")
#  As$LeukemiaType <- relevel(factor(as.character(As$LeukemiaType)), "NoL")
#  A <- subset(A, select = grepl("^ENS", colnames(A)))
#  A <- sapply(A, as.integer)
#  A <- t(A)
#  
#  DEG <- DESeqDataSetFromMatrix(countData = A, colData = As, design = ~ LeukemiaType)
#  DEG <- DESeq(DEG, fitType = "mean", parallel = T, BPPARAM = MulticoreParam(20))
#  
#  n <- levels(As$LeukemiaType)[-1]
#  leukemiasGO <- lapply(setNames(n, n), function (n) {
#    x <- results(DEG, c("LeukemiaType", n, "NoL"), "LeukemiaType", alpha = .05)
#    x <- as.data.frame(x)
#    x$padj[is.na(x$padj)] <- 1
#    g <- new("topGOdata", ontology = "BP", description = 'Leukemia',
#             allGenes = setNames(x$padj, rownames(x)),
#             mapping = "org.Hs.eg.db",
#             geneSelectionFun = function (allScore) { allScore <= 0.05 },
#             annotationFun = annFUN.org, ID = "Ensembl")
#    t <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test") # test definition
#    s <- getSigGroups(g, t) # run F-test
#    r <- GenTable(g, pvalCutOff = s, topNodes = length(g@graph@nodes)) # return data.frame
#    r$pvalCutOff <- as.numeric(str_replace_all(r$pvalCutOff, "[^0-9e\\-\\.]*", ""))
#    r$LogEnrich <- log2(r$Significant / r$Expected)
#    ga <- genesInTerm(g) # GenesAnnotated | list of genes per go-terms
#    ga <- ga[r$GO.ID] # eliminate missing terms
#    names(ga) <- NULL
#    r$GenesAnnotated <- ga
#    xs <- x[,c("padj", "log2FoldChange")] # significant stats subset
#    xs <- subset(xs, padj < 0.05)
#    r$GenesSignificant <- lapply(r$GenesAnnotated, intersect, rownames(xs)) # extract genes
#    ei.rows <- mclapply(r$GenesSignificant, function (y) {
#      if (length(y)) as.list(xs[y,,drop=FALSE])
#      else as.list(rep(NA_real_, length(xs)))
#    }, mc.cores = 10)
#    ei <- mclapply(names(xs), function(z) {
#      lapply(ei.rows, "[[", z)
#    }, mc.cores = 10)
#    ei <- structure(ei, names = names(xs), row.names = seq(nrow(r)), class = "data.frame")
#    row.names(ei) <- NULL
#    r <- data.frame(r, ei, stringsAsFactors = FALSE, check.names = FALSE)
#    return(r)
#  })
#  
#  # lapply(GO, function(x) {list(all(r$Annotated == sapply(r$GenesAnnotated, length)),
#  #                              all(r$Significant == sapply(r$GenesSignificant, length)))})
#  
#  leukemiasGO <- lapply(leukemiasGO, function(x) subset(x, pvalCutOff < 0.1))
#  save(leukemiasGO, file = "data/leukemiasGO.rdata")

## ---- eval=TRUE, echo=FALSE----------------------------------------------
library(devtools)

## ---- eval=FALSE---------------------------------------------------------
#  install_github('dieterich-lab/CellPlot', build_vignettes = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  # load bioconductor functions
#  source("https://bioconductor.org/biocLite.R")
#  
#  # install packages
#  biocLite("BiocParallel")
#  biocLite("annotate")
#  biocLite("DESeq2")
#  biocLite("topGO")
#  biocLite("leukemiasEset")

## ---- eval = FALSE-------------------------------------------------------
#  # install packages
#  install.packages(c("stringr","dplyr","knitr","rmarkdown","devtools"))

## ------------------------------------------------------------------------
# session info from build machine
sessionInfo()

