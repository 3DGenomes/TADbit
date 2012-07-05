batch_tadbit <- function(directory=getwd(), sep='_',
      read_options=list(), verbose=TRUE, ...) {
# Call 'tadbit' in batch on a directory of data files.
# Arguments:
#   directory:    the directory where to find data files.
#   sep:          the string separator for unit names.
#   read_options: arguments passed to 'read.delim'.
#   ...:          further arguments passed to 'tadbit'.
# Return a list of breakpoints per unit. The list names are the
# names of the units found in the data.


   # Here-function for the 'tryCatch' statement (see below).
   tadbit_on_files <- function(fnames) {
      data_list <- list()
      for (fname in fnames) {
         read_options[["file"]] <- fname
         x <- do.call(read.delim, read_options)
         data_list[[fname]] <- as.matrix(x)
      }
      return(tadbit(x=data_list, verbose=verbose, ...))
   }


   # Get file names and unit (chromosome) names.
   fname <- dir(directory, recursive=TRUE, full.names=TRUE)
   unit <- sapply(strsplit(basename(fname), sep, fixe=TRUE), "[", 1)
   fname_list <- tapply(X=fname, INDEX=unit, c)

   # Cycle through units (chromosomes).
   break_list <- list()
   for (unit in names(fname_list)) {
      if (verbose) {
         cat(paste("processing", unit, "\n"))
      }
      # Embed the following in a 'try' for batch robustness.
      tryCatch(
         break_list[[unit]] <- tadbit_on_files(fname_list[[unit]]),
         error = function(e) {
            cat(paste("error with", unit, "(skipping)\n"))
         }
      )
   }

   return (break_list)

}
