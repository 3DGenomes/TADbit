dyn.load("tadbit_R.so")

do_tadbit = function (chr_list) {
   dirs = list.dirs("~/Dropbox/Fcs\ Guillaume/Chromosomes\ 100\ kb/")[-1]
   results = list()
   for (chr_num in chr_list) {
      cat(paste(chr_num), "\n")
      file_name_start = paste("Chr", chr_num, "_", sep="")
      mat_list = list()
      for (directory in dirs) {
         # There should be only one.
         fname = dir(directory, pattern=file_name_start, full.names=TRUE)
         x = read.delim(fname)
         x = as.matrix(x[,-1]) + double(1)
         # Set NAs on rows/columns that have 0 on the diagonal.
         n <- dim(x)[1]
         x <- x[1:(n-2),1:(n-2)]
         to_na = diag(x) < 2
         x[to_na,] = NA
         x[,to_na] = NA
         mat_list[[fname]] = x
      }
      start_time = Sys.time()
      breaks = .Call("tadbit_R_call", mat_list, as.double(50),
            as.integer(10), TRUE)
      print (Sys.time() - start_time)
      results[[chr_num]] = breaks

      # Write to file.
      break_pos = which(breaks == 1)
      start_pos = colnames(x)[c(1, break_pos - 1)]
      end_pos = colnames(x)[c(break_pos, nrow(x))]
      cat(file="out.txt", paste("Chr", chr_num, "\n"), append=TRUE)
      write.table(file="out.txt", sep="\t", quote=FALSE,
         row.names=FALSE, col.names=FALSE,
         data.frame(start=start_pos, end=end_pos),
         append=TRUE)
   }

   return (results)

}

cat("sourced\n")
l = do_tadbit("21")
