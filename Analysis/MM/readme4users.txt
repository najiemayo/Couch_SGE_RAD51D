## to run 
## in R studio console

setwd("where/you/have/downloaded/the/codes/")
## create result folder under the code folder 
if (!dir.exists("results")) {
  dir.create("results")
}
source("./runMM_RAD51D_paper.R")