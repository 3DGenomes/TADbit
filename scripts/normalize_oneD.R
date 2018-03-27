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
cg  = read.csv(args[4], head=FALSE,sep=",")

info = data.frame(t(tot),t(map),t(res),t(cg))
colnames(info) = c("tot","map","res","cg")

info_oned <- oned(info)

write.table(info_oned, file = args[5],row.names = FALSE,col.names = FALSE,sep = ',')

# TODO: use directly the function:
# Compute oned correction
# #'
# #' This function takes a \code{data.frame} with bin information and returns a vector of  biases to correct for
# #' @importFrom mgcv gam
# #' @importFrom mgcv negbin
# #' @importFrom mgcv nb
# #' @param dat A \code{data.frame} with one bin per row containing the total number of contacts and the potential biases as columns.
# #' @param form A \code{formula} describing the  total number of contacts on LHS and the smoothed biases on the RHS (should be compatible with \code{\link{mgcv::gam}})
# #' @return A vector of length \code{nrow(dat)} with the biases to correct for.
# #' @details Please note that the biases returned are the squarerooted so one can directly apply \code{\link{correct_mat_from_b}}
# #' @export
# #' @examples
# #' plot(0)
#
# oned <- function(dat, form = tot ~ s(map) + s(cg) + s(res)){
#
#     fit <- mgcv::gam(as.formula(form), data = dat, family = mgcv::nb())
#
#     out <- predict(fit, newdata = dat, type = "response")
#
#     i_na <- is.na(dat[, all.vars(form)[1]])
#
#     out[which(i_na)] <- NA
#
#     as.numeric(sqrt(out / mean(out, na.rm = T)))
#
# }
#}
