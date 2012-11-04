tadbit <- function(x, n_CPU="auto", verbose=TRUE, max_tad_size="auto",
   heuristic=TRUE) {

   # Validate input type or stop with meaningful error message.
   if (!is.list(x)) {
      if(!is.matrix(x)) {
         stop("'x' must be a matrix or a list of matrices")
      }
      else {
         x <- list(x)
      }
   }
   if (!n_CPU=="auto" && !is.numeric(n_CPU)) {
      stop("'n_CPU' must be \"auto\" or a number")
   }
   if (!max_tad_size == "auto" && !is.numeric(max_tad_size)) {
      stop("'max_tad_size' must be \"auto\" or a number")
   }
   if (!all(sapply(x, is.matrix))) {
      stop("all the elements of 'x' must be matrices")
   }
   ref_dim = dim(x[[1]])
   if (diff(ref_dim)) {
      stop("all the matrices in 'x' must be square")
   }
   for (this_dim in lapply(x, dim)) {
      if (any(this_dim != ref_dim)) {
         stop("all the matrices in 'x' must have same dimensions")
      }
   }
   # Make sure values in 'x' are double.
   # NOTE: calling the C function 'FLOAT' on integers passed from R
   # is not reliable (values differ a bit) and can cause segmentation
   # faults.
   x <- lapply(x, function(y) { y + double(1) })

   # Assign automatic variables and coerce to proper type.
   n_CPU <- as.integer(ifelse(n_CPU == "auto", 0, n_CPU))
   verbose <- as.logical(verbose)
   max_tad_size <- as.integer(ifelse(max_tad_size == "auto",
      ref_dim[1], max_tad_size))
   print(max_tad_size);
   heuristic <- as.integer(heuristic)

   tadbit_c_out <- (.Call("tadbit_R_call", x, n_CPU, verbose,
      max_tad_size, heuristic))

   opt_nbreaks <- tadbit_c_out[[1]]
   llik_mat <- tadbit_c_out[[2]]

   position <- which(tadbit_c_out[[4]][,opt_nbreaks] == 1)
   score <- tadbit_c_out[[5]][tadbit_c_out[[5]] > 0]

   return (list(position=position, score=score))

}
