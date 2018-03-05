#!/usr/bin/env Rscript

library("dryhic")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  stop("At least four arguments must be supplied (tot,map,res,cg).n", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  args[5] = "oneD_biases.csv"
}

tot = read.csv(args[1], head=FALSE,sep=",")
map = read.csv(args[2], head=FALSE,sep=",")
res = read.csv(args[3], head=FALSE,sep=",")
cg = read.csv(args[4], head=FALSE,sep=",")

info = data.frame(t(tot),t(map),t(res),t(cg))
colnames(info) = c("tot","map","res","cg")

info_oned <- oned(info)

write.table(info_oned, file = args[5],row.names = FALSE,col.names = FALSE,sep = ',')