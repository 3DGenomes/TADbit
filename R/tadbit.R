tadbit <- function(x, max_size="auto", n_CPU="auto", verbose=TRUE) {

   # Validate input type or stop with meaningful error message.
   if (!is.list(x)) {
      if(!is.matrix(x)) {
         stop("'x' must be a matrix or a list of matrices")
      }
      else {
         x <- list(x)
      }
   }
   if (!max_size == "auto" && !is.numeric(max_size)) {
      stop("'max_size' must be \"auto\" or a number")
   }
   if (!n_CPU=="auto" && !is.numeric(n_CPU)) {
      stop("'n_CPU' must be \"auto\" or a number")
   }
   if (!all(sapply(x, is.matrix))) {
      stop("all the elements of 'x' must be matrices")
   }
   ref_dim = dim(x[[1]])
   if (diff(ref_dim)) {
      stop("all the matrices in 'x' must be square")
   }
   for (this_dim in lapply(x, dim)) {
      if (this_dim != ref_dim) {
         stop("all the matrices in 'x' must have same dimensions")
      }
   }

   # Assign automatic variables and coerce to proper type.
   max_size <- as.double(ifelse (max_size == "auto", .1, max_size))
   n_CPU <- as.integer(ifelse(n_CPU == "auto", 0, n_CPU))
   verbose <- as.logical(verbose)

   ctdabit <- .Call("tadbit_R_call", x, max_size, n_CPU, verbose)

   return (which(ctadbit == 1))

}
