# We have a function from R^(n x d)
# n is the first elements (n=3 for a triangle)
# d is the dimension of the space (d=2 for R^2)
# For example, an element of R^(n x d) can be (with n=3 and d=2)
# 0 0 
# 1 0
# 1 1
# This means x1=(0,0), x2=(1,0) and x3=(1,1)
# 
# We want to create all possible inputs of R^(n x d) for a 
# specific grid.
# For example, if the grid is -2:2, we want to map for the previous example:
# (-2:2) x (-2:2) x (-2:2) x (-2:2) x (-2:2) x (-2:2),
# which contains 5^6 = 15625 elements.
#
# This can be huge for high resolution, so we first define a shape where some
# elements are already filled. For example, the shape can be:
# NA 0
# NA 0
# NA NA
# This means that x1=(NA,0), x2=(NA,0) and x3=(NA,NA), and we only want
# to fill the 4 NA dimensions with the grid.
#
# The following function do this, with:
# * shape a matrix containing some NA
# * seq_grid the grid for each NA dimension to fill.
input_array = function(shape, seq_grid, limit_size_grid = 100000) {
  colnames(shape) = paste("dim", 1:ncol(shape), sep="")
  rownames(shape) = paste("x", 1:nrow(shape), sep="")
  
  # dim_to_grid corresponds to the NA index to be replaced by seq_grid
  dim_to_grid = which(is.na(shape), arr.ind = TRUE)
  nrow = nrow(dim_to_grid) # additional dimensions to add to the shape
  
  # Test is the shape has no missing elements
  if(nrow == 0) {
    return(shape)
  }
  
  # X is the final array. dim1 and dim2 corresponds to the shape,
  # whereas each higher dimension corresponds to one filled NA
  length_seq_grid = length(seq_grid)
  dim_of_grid = rep(length_seq_grid, nrow)
  X = array(data=rep(shape, prod(dim_of_grid)),
            dim=c(dim(shape), dim_of_grid))
  
  # We want to fill NA in X with the grid.
  # We create the grid:
  grid = grid_matrix(seq_grid, nrow, limit_size_grid)
  
  # Then we fill X with this grid
  for(i in 1:nrow) {
    # We want to write the following:
    # X[dim_to_grid[i,1],dim_to_grid[i,2],,] = grid[,i]
    # But it should work for all hyperdimensions:
    list_to_expand = current_list_to_expand(i, shape, length_seq_grid)
    newvals = expand.grid(list_to_expand)
    newvals = as.matrix(newvals)
    X[newvals] = grid[,i]
  }
  
  # dimnames
  colnames(X) = colnames(shape)
  rownames(X) = rownames(shape)
  names(dimnames(X))[1] = "x"
  names(dimnames(X))[2] = "dim"
  for(i in 1:nrow) {
    fixed_i = paste("x", dim_to_grid[i,1], "[dim", dim_to_grid[i,2], "]", 
                    sep ="")
    dinames_i = paste(fixed_i,
                      "=", seq_grid, sep = "")
    dimnames(X)[[i+2]] = dinames_i
    names(dimnames(X))[i+2] = fixed_i
  }

  return(X)
}

# Shape of the grid, which is the cartesian product of seq_grid over
# the nrow dimensions to fill.
grid_matrix = function(seq_grid, nrow, limit_size_grid = 100000) {
  # Test if the computational cost is reasonable
  size_grid = length(seq_grid)^nrow
  if(size_grid > limit_size_grid) {
    stop(paste("Grid too large. There are ", size_grid, 
               " elements in the grid. You may reduce the grid resolution.",
               sep = ""))
  }
  
  list_to_expand = list()
  for(i in 1:nrow) {
    list_to_expand[[i]] = seq_grid
  }
  
  expanded_grid = expand.grid(list_to_expand)
  return(expanded_grid)
}

# Filling of the output array X with the defined grid
# This function is tricky.
# We want to find all index to fill for each NA in shape.
#
# Each row of dim_to_grid is a value to fill. We take row i.
# For this fixed value, we want to replace all NA in dimensions to fill.
# 
# Example:
# shape = NA 0
#         NA 0
#         NA NA
#
# For i=1, we want to fill the NA on the upper right.
# The two first dimensions is for the shape.
# Each 4 other dimension is for each NA.
current_list_to_expand = function(i, shape, length_seq_grid) {
  dim_to_grid = which(is.na(shape), arr.ind = TRUE)
  nrow = nrow(dim_to_grid)
  ncol = 2 # the shape is always a matrix, ncol(dim_to_grid)
  dim_to_grid_i = dim_to_grid[i,] # the fixed row and column in the shape
  
  list_to_expand = list()
  # The two first dimensions is the shape position
  for(j in 1:ncol) {
    list_to_expand[[j]] = dim_to_grid_i[j]
  }
  # The other dimensions are NA to be replaced by the grid
  for(k in 1:nrow) {
    list_to_expand[[k+ncol]] = 1:length_seq_grid
  }
  return(list_to_expand)
}