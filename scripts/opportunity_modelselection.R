library("callr")
library("Rcpp")
library("RcppArmadillo")
library("tidyverse")
library("Matrix")
library("reshape2")
sourceCpp("scripts/NMF2oppMATSpeed.cpp")

args<-commandArgs(trailingOnly = TRUE) 


type <- args[1]
kmer  <- args[2]
sig <- as.numeric(args[3])
path_to_input <- args[4]
path_to_output <- args[5]

data <- readRDS(path_to_input)

count_matrix <- data[["count_matrix"]]
opportunity_matrix <- data[["background_matrix"]]

init = 3 # number of initialization (we recommend to keep it fixed to 3)
siter = 10000 # number of iteration for the intializations # previously 5000
miter = 100000  #previously 40000

model <-  nmfprmmat(as.matrix(count_matrix), noSignatures = sig, opportunities = as.matrix(opportunity_matrix) ,
                          initial = init, smallIter = siter, maxiter = miter, parametrize = F)

#model_name <- paste(c(paste(c("model", kmer, type), collapse = "_"),"rds"), collapse = ".")

saveRDS(model, file = path_to_output)
