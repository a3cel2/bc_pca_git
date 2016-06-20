library(devtools)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

go_map = '../data/funcassociate_go_associations.txt'
devtools::load_all('../packages/bcPcaAnalysis')
devtools::document('../packages/bcPcaAnalysis')
