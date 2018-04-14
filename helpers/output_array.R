mapX_func = function(X, map_func = map) {
  # apply function for each matrix of the grid
  mapX = apply(X, 3:length(dim(X)), map_func)
  # reshape to have same dim as X
  mapX = array(mapX, dim = dim(X))
  dimnames(mapX) = dimnames(X)
  return(mapX)
}