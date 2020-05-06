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

if(!require(futile.logger)){
  install.packages("futile.logger", dependencies = TRUE)
  library(futile.logger)
}
if(!require(futile.logger)){
  install.packages("xlsx", dependencies = TRUE)
  library(xlsx)
}


quit(save="no")