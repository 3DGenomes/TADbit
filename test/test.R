
library(tadbit)

## R --no-save --slave < test.R

## TEST TADBIT
x = as.matrix(read.delim('chrT/chrT_A.tsv', row.names=1))
output_fast = tadbit(x)
exp.pos <- c(8, 14, 19, 34, 39, 44, 50, 62, 67, 72, 90)
exp.scr <- c(3, 6, 7, 4, 4, 8, 9, 6, 8, 7, 9)
if ((output_fast$position==exp.pos) && (output_fast$score==exp.scr)) {
  print("Passed tadbit test...")
}else{
  stop("Error unexpected result for tadbit")
}

## TEST BATCH_TADBIT
output = batch_tadbit('chrT/', max_tad_size=4, no_heuristic=TRUE)
exp.pos <- c(5, 14, 19, 34, 39, 47, 53, 62, 67, 72, 90, 95)
exp.scr <- c(3, 6, 8, 8, 8, 5, 4, 5, 5, 7, 7, 7)
if ((output$chrT$position==exp.pos) && (output$chrT$score==exp.scr)) {
  print("Passed tadbit test...")
}else{
  stop("Error unexpected result for tadbit")
}

print("all tests ok!")
