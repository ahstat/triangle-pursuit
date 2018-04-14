##
# Translation
##
# Translation of vector t
formula_translation= function(shape, t) {
  return(shape + t)
}

# Homothety to (0,0) with factor t
formula_homothety= function(shape, t) {
  return(shape * t)
}

# Rotation to (0,0) with angle t
formula_rotation = function(shape, t) {
  if(ncol(shape) != 2) {
    stop("Rotation coded only for 2 dimensional polygons")
  }
  complex_out = (shape[,1]+1i*shape[,2])*exp(1i*t)
  output = matrix(NA, nrow = nrow(shape), ncol = ncol(shape))
  output[,1] = Re(complex_out)
  output[,2] = Im(complex_out)
  return(output)
}

input_operation = function(shape, seq_grid, name_operation) {
  colnames(shape) = paste("dim", 1:ncol(shape), sep="")
  rownames(shape) = paste("x", 1:nrow(shape), sep="")
  
  # X is the final array. dim1 and dim2 corresponds to the shape,
  # whereas each higher dimension corresponds to one filled NA
  length_seq_grid = length(seq_grid)
  X = array(data=rep(NA, length(shape)*length_seq_grid),
            dim=c(dim(shape), length_seq_grid))
  
  formula_operation = NA
  if(name_operation == "translation") {
    formula_operation = formula_translation
  } else if(name_operation == "homothety") {
    formula_operation = formula_homothety
  } else if(name_operation == "rotation") {
    formula_operation = formula_rotation
  }
  
  for(j in 1:length_seq_grid) {
    X[,,j] = formula_operation(shape, seq_grid[j])
  }
  
  # dimnames
  colnames(X) = colnames(shape)
  rownames(X) = rownames(shape)
  dimnames(X)[[3]] = paste(name_operation, seq_grid)
  names(dimnames(X))[1] = "x"
  names(dimnames(X))[2] = "dim"
  names(dimnames(X))[3] = name_operation

  return(X)
}