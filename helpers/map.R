# From a vector of size k, deduce a mapped vector which is also of size k
map = function(x, ...) {
  
  if(is.null(dim(x))) {
    # convert to a matrix if we get a vector.
    # (in dim 1; or in dim 2 with a vector of complex numbers)
    warning("The input should be a column matrix, not a vector.")
    x = as.matrix(x)
  }
  len = dim(x)[1]

  seq = get_recursive(x, ...)
  
  if(any(is.na(seq))) { # if we cannot compute the whole sequence
    map_x = matrix(NA, nrow = nrow(x), ncol = ncol(x))
  } else {
    last_elements = tail(seq, len)
    
    # ordering with corresponding vertex
    # For example, if the initial sequence is x[1], x[2], x[3],
    # then, x[4] is in the direction x[1], x[5] in the direction x[2], etc.
    # Finally, the converged x[n], x[n+1], x[n+2] are in directions
    # x[a], x[b], x[c], and we find {p,q,r}={n,n+1,n+2} such that:
    # x[p] is in direction x[1],
    # x[q] is in direction x[2],
    # x[r] is in direction x[3].
    # The output map is: (x[p], x[q], x[r])
    last_indexes = tail(1:nrow(seq), len)
    modulo = last_indexes %mod% len
    map_x = last_elements[order(modulo),,drop = FALSE]
  }
  
  return(map_x)
}

# Function which output a custom map function
map_custom = function(nwalk = 1000, 
                      type_norm = euc_norm, 
                      epsilon = 1e-8,
                      rule = rule_shift) {
  func = function(x) {
    map(x, nwalk = nwalk, type_norm = type_norm, 
        epsilon = epsilon, rule = rule)
  }
  return(func)
}

# Plotting function for map
# circle for inputs, square for outputs
# Color order: black, red, green, blue, cyan... 
map_plot = function(x, ...) {
  y = map(x, ...)
  
  xlab = "x"
  ylab = "y"
  
  if(is.complex(x)) {
    xlim = range(Re(c(x,y)))
    ylim = range(Im(c(x,y)))
  } else if(ncol(x)==1) {
    xlim = range(c(x,y))
    
    # dimension y only to distinguish x and map(x)
    ylim = c(0, 1)
    ylab = "0 for input, 1 for output"
    x = cbind(x, rep(0, length(x)))
    y = cbind(y, rep(1, length(y)))
  } else {
    xlim = range(rbind(x,y)[,1]) 
    ylim = range(rbind(x,y)[,2])
  }

  col =1:nrow(x)
  plot(x, col = col, 
       xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab)
  lines(y, col = col, type = "p", pch = 0)
}