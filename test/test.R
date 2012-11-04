


library(tadbit)
# 'x' is a square matrix with raw counts. It can also be a list of squre matrices with identical dimensions.
x = as.matrix(read.delim('chrT/chrT_A.tsv', row.names=1))
# Plot a heatmap of the hiC map.
# Call tadbit
# output_slow = tadbit(x)
output_fast = tadbit(x, max_tad_size=4, heuristic=TRUE)
print(output_fast)

image(log2(x))
for (l in output_fast[1]){
  abline(h=l/100)
  abline(v=l/100)
}

# You can also call a batch tadbit on all 4 files at the same time
batch_output = batch_tadbit('chrT/', max_tad_size=4,  n_CPU=1, heuristic=TRUE,
  read_options=list(row.names=1))
print (batch_output)

image(log2(x))
for (l in batch_output$chrT[1]){
  print (l)
  abline(h=l/100)
  abline(v=l/100)
}
