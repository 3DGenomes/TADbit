#!/usr/bin/env Rscript

# example usages:
# $ Rscript --vanilla normalize_oneD.R input.csv
# $ Rscript --vanilla normalize_oneD.R input.csv output.tsv  0.1
# $ Rscript --vanilla normalize_oneD.R input.csv output.tsv
# $ Rscript --vanilla normalize_oneD.R input.csv 0.1
# $ Rscript --vanilla normalize_oneD.R input.csv output.tsv "tot ~ s(map) + s(cg)" 0.1

library("dryhic")

# parse arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("At least two arguments must be supplied (tot,map,res,cg).n", call.=FALSE)
}

output <- "oneD_biases.csv"
p_fit <- NA
form <- "tot ~ s(map) + s(cg) + s(res)"

if (length(args) > 1){
    for (i in seq(2, length(args), 1)) {
        if (grepl("~", args[i])){
            form <- args[i]
        }else if(!suppressWarnings(is.na(as.numeric(args[i])))){
            p_fit <- as.numeric(args[i])
        }else{
            output <- args[i]
        }
    }
}

form <- as.formula(form)

# get data
info <- read.csv(args[1],sep=",")

# run oneD
info_oned <- oned(info, form=form, p_fit=p_fit)

# write results
write.table(info_oned, file=output, row.names=FALSE,
            col.names=FALSE, sep = ',')
