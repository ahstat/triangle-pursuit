##
# Input x for get_recursive and related functions.
#
# It should be a matrix of size n x d.
# The sequence must begin with at least n=2 elements x_1, x_2,
# Each element being a row of the matrix x.
# Then, the number of column d is the state space of each element (R^d)
##

# Compute sequence x[1,2,3,...] until convergence, which is recursively
#  defined according to the rule.
# Convergence is defined through epsilon, and we throw an error if there is
#  no convergence after nwalk iterations.
get_recursive = function(x, nwalk = 1000, 
                         type_norm = euc_norm, epsilon = 1e-8,
                         rule = rule_shift) {
  if(is.null(dim(x))) {
    # convert to a matrix if we get a vector.
    # (in dim 1; or in dim 2 with a vector of complex numbers)
    warning("The input should be a column matrix, not a vector.")
    x = as.matrix(x)
  }

  len = dim(x)[1]
  ncol = dim(x)[2]
  
  if(len<=1) {
    stop("The initial sequence should have at least 2 elements.
         Maybe you need to transpose the matrix?")
  }
  
  # output matrix, x[i,] being the term i of the sequence
  x = rbind(x, matrix(rep(NA, ncol*(nwalk-len)), ncol = ncol))
  remain_steps = NA
  for(i in (len+1):nwalk) {
    x[i,] = rule(x, i, len, type_norm)

    # If we want to go to direction 'a' from 'a', then the result
    #   is NaN and we cannot continue iterations.
    if(any(is.nan(x[i,]))) {
      break()
    }
    
    # Check if there is convergence at this step
    if(is_need_stop(x[i-len,], x[i,], type_norm, epsilon)) {
      #break()
      # need to continue, e.g. for (0,1,2,-1) initially
      if(is.na(remain_steps)) {
        remain_steps = len+1 # can be optimized... len? len-1? 
      } else if(remain_steps > 0){
        remain_steps = remain_steps - 1
      } else {
        break()
      }
    } else {
      remain_steps = NA
    }
  }
  
  # If the loop ended without any break, it means that convergence
  #   has not come.
  if(i == nwalk) {
    print("Current x:")
    print(x[1:20,,drop = FALSE])
    stop(paste("nwalk is too small, need to increase it...
         Or not", len, "adherent points."))
  }
  
  # Only return the sequence before convergence
  x = x[1:i,,drop = FALSE]
  return(x)
}

# Test to know whether we stop the for loop now
is_need_stop = function(x, y, type_norm, epsilon) {
  return(ifelse(type_norm(x-y) < epsilon, TRUE, FALSE))
}

# Plotting function for get_recursive, just a test function
get_recursive_plot = function(x, only_tails = FALSE, ...) {
  y = get_recursive(x, ...)
  
  if(only_tails) {
    len = dim(x)[1]
    y = tail(y, len+1)
  }
  
  plot(y, type = 'o', asp = 1, xlab = "x", ylab = "y")
}