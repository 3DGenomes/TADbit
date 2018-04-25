#!/usr/bin/env Rscript

# example usages:
# # only input file
# $ Rscript --vanilla normalize_oneD.R input.csv
# # input file and using 10% of the data
# $ Rscript --vanilla normalize_oneD.R input.csv output.tsv  0.1
# # input file and output file
# $ Rscript --vanilla normalize_oneD.R input.csv output.tsv
# # input file, using 10% of the data, and setting the seed number to 2 (input seed should always be > 1)
# $ Rscript --vanilla normalize_oneD.R input.csv 0.1 2
# # input file, output file, formula and 10% of the input data used
# $ Rscript --vanilla normalize_oneD.R input.csv output.tsv "tot ~ s(map) + s(cg)" 0.1

library("dryhic")

# parse arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("At least two arguments must be supplied (tot,map,res,cg).n", call.=FALSE)
}

output <- "oneD_biases.csv"
p_fit <- NA
seed <- 1
form <- "tot ~ s(map) + s(cg) + s(res)"

if (length(args) > 1){
    for (i in seq(2, length(args), 1)) {
        if (grepl("~", args[i])){
            form <- args[i]
        }else if(!suppressWarnings(is.na(as.numeric(args[i])))){
            if(as.numeric(args[i]) <= 1){
                p_fit <- as.numeric(args[i])
            }else{
                seed <- as.numeric(args[i])
            }

        }else{
            output <- args[i]
        }
    }
}

form <- as.formula(form)

# get data
info <- read.csv(args[1],sep=",")

# initialize seed
set.seed(seed)

# run oneD
info_oned <- oned(info, form=form, p_fit=p_fit)

# write results
write.table(info_oned, file=output, row.names=FALSE,
            col.names=FALSE, sep = ',')
