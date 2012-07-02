batch_tadbit <- function(directory=getwd(), sep='_',
      read_options=list(), ...) {
# Call |tadbit| in batch on a directory of data files.
# Arguments:
#   directory:    the directory where to find data files.
#   sep:          the string separator for unit names.
#   read_options: arguments passed to |read.delim|.
#   ...:          further arguments passed to |tadbit|.
# Return a list of breakpoints per unit. The list names are the
# names of the units found in the data.

   # Get file names and unit (chromosome) names.
   fname <- dir(directory, recursive=TRUE, full.names=TRUE)
   unit <- sapply(strsplit(basename(fname), sep), "[", 1)
   fname_list <- tapply(X=fname, INDEX=unit, c)

   # Cycle through units (chromosomes).
   break_list <- list()
   for (unit_name in names(fname_list)) {
      data_list <- list()
      for (fname in fname_list[[unit_name]]) {
         read_options[["file"]] <- fname
         x <- do.call(read.delim, read_options)
         data_list[[fname]] <- as.matrix(x)
      }
      break_list[[unit_name]] <- tadbit(x=data_list, ...)
   }

   return (break_list)
}
