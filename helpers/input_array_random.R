input_random = function(shape, nb_set, val_min=-1, val_max=1) {
  colnames(shape) = paste("dim", 1:ncol(shape), sep="")
  rownames(shape) = paste("x", 1:nrow(shape), sep="")
  
  # X is the final array. dim1 and dim2 corresponds to the shape
  X = array(data=rep(shape, length(shape)*nb_set),
            dim=c(dim(shape), nb_set))
  
  to_replace = which(is.na(X))
  X[to_replace] = runif(length(to_replace), val_min, val_max)
  
  # dimnames
  colnames(X) = colnames(shape)
  rownames(X) = rownames(shape)
  dimnames(X)[[3]] = paste("random", 1:nb_set)
  names(dimnames(X))[1] = "x"
  names(dimnames(X))[2] = "dim"
  names(dimnames(X))[3] = paste("runif(", val_min, ", ", val_max, ")", sep = "")
  
  return(X)
}