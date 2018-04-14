###
# "Triangle pursuit" project.
# (2017/05/20)
#
# The initial problem is as follows: given 3 points in R^2, called
# x1, x2, x3, we define the recurrent sequence with:
#     x_i = x_{i-1} + N(x_{i-1} - x_{i-3}),
# where N(.) is the unit normalization defined with, for x in R^2:
#     N(x) = x / ||x||,
# with ||.|| the Euclidian norm.
#
# The general problem is introduced here and in the following paragraphs.
# Given 'len' initial points of R^d (with 'len'>1 and d>0), we define 
# the recurrent sequence with:
#     x_i = rule(x_{1:(i-1)}),
# with a rule with the same shape as for the initial problem.
#
# The sequence is coded with a matrix x, the columns representing the dimension
# and the rows representing the elements x_i. For example, with 3 points in
# R^2, the 3 initial points are coded as as 3x2 matrix:
#     x_1[1] x_1[2]
#     x_2[1] x_2[2]
#     x_3[1] x_3[2]
##
rm(list = ls())
setwd("~/Documents/GitHub/triangle-pursuit/")
verbose = TRUE

#########
# Rules #
#########
##
# First we define:
# * a norm, for example ||.||,
# * a rule, defining how to compute x_i from x_{1:(i-1)}
##
source("helpers/rules.R")

##
# Here is the available norms.
# For vect an element of R^d (with d a positive integer), we have:
# * Euclidian norm: euc_norm(vect)
# * Max norm: max_norm(vect)
# * One norm: one_norm(vect)
##

##
# Assuming the sequence has 'len' initial points, 'type_norm' the selected
# norm, x representing the i-1 first elements of the sequence, we define
# the i-th element with:
#     rule(x, i, len, type_norm)
#
# Two rules are available:
# * rule_shift (default): x_i = x_{i-1} + N(x_{i-1} - x_{i-len}),
# * rule_origin: x_i = x_{i-1} + N(x_{i-1} - x_{i modulo len}),
# where i modulo len is defined to have an element in x_1, ... x_{len}.
##

############
# Sequence #
############
##
# Then we compute the sequence x_{i} until convergence to its adherent points.
# We assume there are 'len' adherent points, which is however not true with all
# defined norms and rules.
##
source("helpers/get_recursive.R")

##
# The function to compute the sequence is of the form:
#     get_recursive(x, nwalk = 1000, type_norm = euc_norm, 
#                   epsilon = 1e-8, rule = rule_shift) 
# * x is the initial matrix of size 'len' x m,
# * nwalk is the maximum number of steps computed before throwing an error,
# * type_norm is the type of the norm,
# * epsilon is the precision, i.e. we stop when |x_{i}-x_{i-len}|<epsilon is
# reached for an entire sequence of length len,
# * rule is the selected rule.
##

# Simple tests for the get_recursive function
if(verbose) {
  ##
  # With the initial triangle x1 = (0,0), x2 = (1,0), x3 = (1,1) in 2D
  ##
  x = rbind(c(0,0), c(1,0), c(1,1))
  # x = as.matrix(c(0,1,1+1i)) # Works with complex numbers, but experimental
  
  # Trajectory
  get_recursive_plot(x, rule = rule_shift)
  
  # Final triangle
  get_recursive_plot(x, rule = rule_shift, only_tails = TRUE)
  
  # Final triangle, changing the rule to rule_origin
  # We obtain a slightly different triangle
  get_recursive_plot(x, rule = rule_origin, only_tails = TRUE)
  
  ##
  # With the initial quadrilater x1=(0,0), x2=(1,0), x3=(1,1), x4=(0,2)
  ##
  x = rbind(c(0,0), c(1,0), c(1,1), c(0,2))
  # x = as.matrix(c(0,1,1+1i,2i))
  
  # Trajectory
  get_recursive_plot(x, rule = rule_shift)
  
  # Final quadrilater
  get_recursive_plot(x, rule = rule_shift, only_tails = TRUE)
  
  ##
  # With the initial triangle x1=(0,0,0), x2=(1,0,0), x3=(1,1,1)
  ##
  x = rbind(c(0,0,0), c(1,0,0), c(1,1,1))
  y = get_recursive(x)
  print("===Trajectory from x1=(0,0,0), x2=(1,0,0), x3=(1,1,1)===")
  print(y)
}

# Tricky tests for the get_recursive function
if(verbose) {
  # Correct output for this initialization:
  print(get_recursive(t(t(c(0,1,2,-1)))))
  
  # The following throws an error
  # This is beause there are 6 adherent points instead of 3
  # (with the rule_shift rule, I have not experienced this problem).
  x = matrix(c(0,0,1,-0.1,2,1), nrow = 3, ncol = 2, byrow = TRUE)
  #get_recursive(x, rule = rule_origin) # error
  ##in the debug mode, this can highlight there are 6 adherent pts:
  #plot(x[,1]) 
  #plot(x[,2])
}

###########
# Mapping #
###########
##
# From x_1, ..., x_len, we obtained the sequence until convergence.
# In the map function, we deduce the 'len' adherent points.
# This means that: map(x_1, ..., x_len) is also a vector of 'len' points called
# y_1, ..., y_len.
#
# y_i is linked with x_i, which means that y_i = lim_n x_{len*n + i}.
##
source("helpers/map.R")

##
# The function for mapping is of the form:
#     map(x, ...),
# where: 
# * x is the initial matrix of size 'len' x m,
# * ... can be any option of get_recursive.
#
# Custom maps can be created with the map_custom function:
#     map_custom(nwalk = 1000, type_norm = euc_norm,
#                epsilon = 1e-8, rule = rule_shift)
##

# Simple tests for the map function
if(verbose) {
  ##
  # With the initial triangle x1 = (0,0), x2 = (1,0), x3 = (1,1) in 2D
  ##
  x = rbind(c(0,0), c(1,0), c(1,1))
  # x = as.matrix(c(0,1,1+1i)) # Working with complex numbers but experimental
  
  # Trajectory with default rule
  # Points correspond to the initial sequence, and squares to the output 
  # mapping.
  print("===Trajectory from x1=(0,0), x2=(1,0), x3=(1,1)===")
  print(map(x))
  map_plot(x)
  
  # Trajectory with the rule_origin rule
  map_plot(x, rule = rule_origin)
  
  # In complex, need to provide the vector as a column matrix
  # x = c(0,1,1+1i)
  # map_plot(x) # not working, because need a column matrix
  
  ##
  # With the two points x1=(0), x2=(0.5) in dimension 1
  ##
  # On the plot, the input is located on the line y = 0 (with points), 
  # and the output on the line y = 1 (with squares).
  x = rbind(c(0), c(0.5))
  y = map_plot(x, rule = rule_origin)
  
  ##
  # With the initial quadrilater x1=(0,0), x2=(1,0), x3=(0,1), x4=(1,-1)
  ##
  x = rbind(0, 1, 1i, 1i-1)
  map_plot(x)
}

####################
# Input arrays (1) #
####################
##
# Instead of an initial matrix x containing the first 'len' elements, and
# of size len x m, we define a shape.
# A shape is a matrix of size len x m, but containing some NA values.
# Those NA values will be replaced by a grid of values.
#
# For example, with the shape:
#             0  0
#     shape = 0  1
#             2 NA
#
# There is one value to fill.
# By defining a grid seq_grid, for example: c(-2,-1,0,1,2), we obtain an
# array X defined with length(seq_grid)^(#NA in the shape) matrices:
#         0  0   0  0   0  0   0  0   0  0
#     X = 0  1   0  1   0  1   0  1   0  1
#         2 -2   2 -1   0  0   0  1   0  2
#
# The function to do this is input_array:
#     input_array(shape, seq_grid, limit_size_grid = 100000)
# Here, the function will stop if the size of the grid is beyond 
# limit_size_grid
source("helpers/input_array.R")

# Tests for input_array, with different shape of matrix
if(verbose) {
  ##
  # Tests with triangle in dimension 2
  ##
  # 'Shape' of the shape
  shape = matrix(NA, nrow = 3, ncol = 2, byrow = TRUE)
  
  # Grid for each NA dimension
  seq_grid = seq(from = -2, to = 2, length.out = 11)
  
  # Example 1
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_array(shape, seq_grid)
  
  # Example 2
  shape = matrix(c(0,0,1,0,1,NA), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_array(shape, seq_grid)
  
  # Example 3
  shape = matrix(c(NA,0,1,0,1,NA), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_array(shape, seq_grid)
  
  # Example 4
  shape = matrix(c(NA,0,NA,0,1,NA), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_array(shape, seq_grid)
  
  # Example 5
  shape = matrix(c(NA,NA,NA,NA,1,NA), nrow = 3, ncol = 2, byrow = TRUE)
  # X = input_array(shape, seq_grid) # Error grid too large
  X = input_array(shape, seq_grid, limit_size_grid = 200000)
  
  ##
  # Tests with two points in dimension 2
  ##
  # Example 6
  shape = matrix(c(1,NA,NA,NA), nrow = 2, ncol = 2, byrow = TRUE)
  X = input_array(shape, seq_grid)
  
  # Example 7
  shape = matrix(c(NA,NA,NA,NA), nrow = 2, ncol = 2, byrow = TRUE)
  X = input_array(shape, seq_grid)
}

####################
# Input arrays (2) #
####################
##
# Instead of a grid, we can be interested by applying an operation on a
# fixed matrix x.
# Three operations are coded:
# * Translation,
# * Homothety,
# * Rotation.
#
# Finally, we obtain an array X with 3 directions.
#
# The function to do this is:
#     input_operation(shape, seq_grid, name_operation)
#
# Where name_operation can be:
# * "translation"
# * "homothety"
# * "rotation"
#
# And where the shape is a matrix which do not contain any NA.
##
source("helpers/input_array_operations.R")

# Some examples for input_operation
if(verbose) {
  # Translation of each coordinate by the same constant
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -2, to = 2, length.out = 11)
  X = input_operation(shape, seq_grid, "translation")
  
  # Rotation around 0
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = 0, to = 2*pi, length.out = 11)
  X = input_operation(shape, seq_grid, "rotation")
  
  # Homothety
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -2, to = 2, length.out = 11)
  X = input_operation(shape, seq_grid, "homothety")
}

####################
# Input arrays (3) #
####################
##
# Sometimes, we may want to let randomness fill the NA of the shape.
#
# This is done with the function:
#     input_random(shape, nb_set, val_min=-1, val_max=1)
# With:
# * shape: a matrix containing some NA to be filled,
# * nb_set: number of times the matrix will be filled,
# * val_min: minimum value of the random distribution,
# * val_max: maximum value of the random distribution.
#
# The distribution is uniform. Finally, we have an array with 3 directions.
##
source("helpers/input_array_random.R")

# Example for the function input_random, filling all NA of the shape
# with random values.
if(verbose) {
  shape = matrix(c(0,0,NA,0,NA,NA), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 5, -1, 1)
  
  shape = matrix(c(NA,NA,NA,NA,NA,NA), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 5, 0, 1)
  
  shape = matrix(c(NA,3,1,0,0,0), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 5, 0, 1)
}

#################
# Array mapping #
#################
##
# Now, we have an array X, each element being a matrix x.
# We can apply the mapping on each element of the array.
# This is done with the function:
#     mapX_func(X, map_func = map)
# with:
# * X an array of inputs,
# * map a map created before.
##
source("helpers/output_array.R")

# A simple example of mapping
if(verbose) {
  # Triangle in dimension 2
  shape = matrix(c(0,0,1,0,NA,NA), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -1, to = 1, length.out = 11)
  X = input_array(shape, seq_grid)
  
  # Check one matrix and the related mapping
  X[,,"x3[dim1]=1","x3[dim2]=1"] # X[,,11,11]
  map(X[,,"x3[dim1]=1","x3[dim2]=1"]) # map(X[,,11,11])

  # Mapping for all matrices of X
  mapX = mapX_func(X)
  dim(mapX)
  
  # Retrive the previous mapping from here
  mapX[,,"x3[dim1]=1","x3[dim2]=1"]
}

#########
# Plots #
#########
##
# We want to represent this mapping.
# Here, the initial shape must have only one NA, so the array X will have 3
# directions (x_i, dim and seq_grid).
#
# In that case, there is a unique plot function:
#
#     plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 1,
#               cols = NA, nb_text = NA, nb_polygon = NA,
#               main = NA, cex = 0.8, type = "p", , ...)
#
# with:
# * seq_grid: vector of the grid (or any vector with the same size),
# * mapX: the mapping obtained with mapX_func(X, map_func = map),
# * xaxis: the x axis. NA to have seq_grid has xaxis, else can be 1 or 2,
# * yaxis: the y axis. NA to have seq_grid has xaxis, else can be 1 or 2,
# * cols: select the color for each output sequence (NA takes the default),
# * nb_text: number of indication of seq_grid plotted on the 2D graph,
# * nb_polygon: number of polygon linking the sequences on the 2D graph,
# * main: plotting parameter, title,
# * cex: plotting parameter, size of text if any,
# * type: plotting parameter, type of points,
# * ...: other plotting parameters.
##
source("helpers/plots.R")

######################
# Plots examples (1) # (input_array)
######################
# Many example of plots with input_array and mapX_func
if(verbose) {
  ##
  # Example 1 (triangle, 2 dimensional space, 1 moving value)
  ##
  # 0 0
  # 1 x
  # 2 1
  shape = matrix(c(0,0,1,NA,2,1), nrow = 3, ncol = 2, byrow = TRUE)
  #seq_grid = seq(from = -2, to = 2, length.out = 2001)
  seq_grid = seq(from = 0, to = 1, length.out = 2001)
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  
  main = paste(as.vector(t(cbind("[", shape, "] "))), collapse = "")
  # Plot the first dim of mapX as a function of seq_grid
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 1, main = main)
  # Plot the second dim of mapX as a function of seq_grid
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 2, main = main)
  # Plot the 2 dimensions together
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, main = main,
            nb_polygon = 2, nb_text = 2, asp = 1)
  # The same plot with lines instead of points
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, main = main, 
            type = "l", nb_text = 2, nb_polygon = 2, asp = 1)
  
  # Plot of the initial triangle x_1, x_2, x_3 moving
  # (to show it does not seem so natural to move one coordinate like this,
  # and there is something implying the discontinuity)
  plot_mapX(seq_grid, X, xaxis = 1, yaxis = 2, main = main, nb_polygon = 10)
  
  ##
  # Example 2 (triangle, 2 dimensional space, 1 moving value, rule = origin)
  ##
  # 0 0
  # 1 x
  # 2 1
  # The grid is only in [-2,-0.2] because the number of adherent points
  # switches from 3 to 6 after that (for example for -0.1).
  shape = matrix(c(0,0,1,NA,2,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -2, to = -0.2, length.out = 2001)
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X, map_func = map_custom(rule = rule_origin,
                                            epsilon = 1e-5,
                                            nwalk = 10000))
  
  main = paste(as.vector(t(cbind("[", shape, "] "))), collapse = "")
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 1, main = main)
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 2, main = main)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, main = main)
  
  ##
  # Example 3 (triangle, 2 dimensional space, 1 moving value)
  ##
  # 0 0
  # 1 0
  # 4 x
  var1 = 1 # any positive number
  var2 = 4 # any number
  var3 = NA #2 # any number
  shape = matrix(c(0,0,var1,0,var2,var3), nrow = 3, ncol = 2, byrow = TRUE)
  
  seq_grid = seq(from = -10, to = 10, length.out = 2001)
  
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  
  main = paste(as.vector(t(cbind("[", shape, "] "))), collapse = "")
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, main = main)
  
  ##
  # Example 4 (quadrilateral, 2 dimensional space, 1 moving value) !!
  ##
  # 0 0
  # 1 0
  # 1 2
  # 2 x
  shape = matrix(c(0,0,1,0,1,2,2,NA), nrow = 4, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -2, to = 2, length.out = 2001)
  
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  
  # With this examples, there are 2 phases:
  # Before a certain value, the final polygon is a quadrilater,
  # After this value, it is 2 segments, and the sequence is of the form
  # (a, b, a, c).
  main = paste(as.vector(t(cbind("[", shape, "] "))), collapse = "")
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 1, main = main)
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 2, main = main)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, main = main, 
            cex = 0.5, nb_text = 10, nb_polygon = 5)
  
  ##
  # Example 5 (pentagon, 2 dimensional space, 1 moving value)
  ##
  # just another plot
  shape = matrix(c(0,0,1,0,1,2,3,2,2,NA), nrow = 5, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -2, to = 2, length.out = 2001)
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.5)
  
  ##
  # Example 6 (hexagon, 2 dimensional space, 1 moving value)
  ##
  # just another plot
  shape = matrix(c(0,0,1,0,1,2,3,2,4,1,2,NA), nrow = 6, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -2, to = 2, length.out = 2001)
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.5)
  
  ##
  # Example 7 (pentagon, 3 dimensional space, 1 moving value)
  ##
  shape = matrix(c(0,0,0,1,0,0,1,2,1,3,2,1,2,1,NA), 
                 nrow = 5, ncol = 3, byrow = TRUE)
  seq_grid = seq(from = -2, to = 2, length.out = 2001)
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.5)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 3, cex = 0.5)
  plot_mapX(seq_grid, mapX, xaxis = 2, yaxis = 3, cex = 0.5)
  
  ##
  # Example 8 (quadrilateral, 3 dimensional space, 1 moving value)
  ##
  # just another plot
  shape = matrix(c(0,0,0,1,0,0,1,2,1,2,1,NA), nrow = 4, ncol = 3, byrow = TRUE)
  seq_grid = seq(from = -2, to = 2, length.out = 2001)
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 3)
  plot_mapX(seq_grid, mapX, xaxis = 2, yaxis = 3)
  
  ##
  # Example 9 (triangle, 3 dimensional space, 1 moving value)
  ##
  shape = matrix(c(0,0,0,1,0,0,0,1,NA), nrow = 3, ncol = 3, byrow = TRUE)
  seq_grid = seq(from = -2, to = 2, length.out = 2001)
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 3)
  plot_mapX(seq_grid, mapX, xaxis = 2, yaxis = 3)
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 1)
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 2)
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 3)
  
  ##
  # Example 10 (triangle, 2 dimensional space, 1 moving value)
  ##
  # 0 0
  # 1 0
  # x 2
  shape = matrix(c(0,0,1,0,NA,2), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -20, to = 20, length.out = 201)
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 1)
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 2)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2)
  
  ##
  # Example 11 (bigone, 1 dimensional space, 1 moving value)
  ##
  # 0
  # x
  #
  # This example can be easily solved analytically
  shape = matrix(c(0,NA), nrow = 2, ncol = 1, byrow = TRUE)
  
  seq_grid = seq(from = -20, to = 20, length.out = 201)
  
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 1)
  
  ##
  # Example 12 (quadrilateral, 1 dimensional space, 1 moving value)
  ##
  # 0
  # 1
  # 2
  # x
  #
  # This can also be solved analytically
  shape = matrix(c(0,1,2,NA), nrow = 4, ncol = 1, byrow = TRUE)
  seq_grid = seq(from = -5, to = 5, length.out = 2001)
  
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 1, cex = 0.1)
  # The plot is not so neat, because we may want to distinguish when
  # two elements share the same values.
  
  ##
  # Example 13 (quadrilateral, 1 dimensional space, 1 moving value)
  ##
  # 0
  # -3.2
  # 4.5
  # x
  shape = matrix(c(0,-3.2,4.5, NA), nrow = 4, ncol = 1, byrow = TRUE)
  
  seq_grid = seq(from = -5, to = 5, length.out = 2001)
  
  X = input_array(shape, seq_grid)
  mapX = mapX_func(X)
  
  plot_mapX(seq_grid, mapX, xaxis = NA, yaxis = 1, cex = 0.1)
}

######################
# Plots examples (2) # (input_operation)
######################
# Examples of plots with input_operation and mapX_func
if(verbose) {
  # Translation keeps the transformation
  # (T*map = map*T)
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -2, to = 2, length.out = 501)
  X = input_operation(shape, seq_grid, "translation")
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 2, nb_text = 5, asp = 1)
  
  # Rotation keeps the transformation (for the Euclidian norm!)
  # (R*map = map*R)
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = 0, to = 2*pi, length.out = 501)
  X = input_operation(shape, seq_grid, "rotation")
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 10, nb_text = 5, asp = 1)
  
  # Homothety does not keep the transformation
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -5, to = 5, length.out = 501)
  X = input_operation(shape, seq_grid, "homothety")
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 10, nb_text = 5, asp = 1)
  
  # Homothety with another moving parameter
  # (just for illustration purpose)
  for(t in seq(from = -1, to = 1, length.out = 50)) {
    shape = matrix(c(1,0,-1,t,0,0), nrow = 3, ncol = 2, byrow = TRUE)
    seq_grid = seq(from = -2, to = 2, length.out = 101)
    X = input_operation(shape, seq_grid, "homothety")
    mapX = mapX_func(X)
    plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p",
              main = t)
  }
  
  # Another homothety with another moving parameter
  # (just for illustration purpose)
  for(t in seq(from = -1, to = 0.99, length.out = 50)) {
    shape = matrix(c(0,0,t,1,1,1), nrow = 3, ncol = 2, byrow = TRUE)
    seq_grid = seq(from = -2, to = 2, length.out = 101)
    X = input_operation(shape, seq_grid, "homothety")
    mapX = mapX_func(X)
    plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p",
              main = t)
  }
}

# Homothety of an equilateral triangle (beautiful to plot) !!
# and calculation for the terms x4, x5, x6.
if(verbose) {
  shape = matrix(c(1,0,1/2,sqrt(3)/2,0,0), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -1, to = 1, length.out = 2001)
  X = input_operation(shape, seq_grid, "homothety")
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 2, nb_text = 5, asp = 1)
  
  # This corresponds to the initialization (with t varying):
  # x1 = (t, 0); x2 = ((1/2)*t, (\sqrt{3}/2)*t); x3 = (0, 0),
  # As seen here for t varying from -1 to 1:
  plot_mapX(seq_grid, X, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 20, asp = 1)
  
  # There is a symmetry, so we can only consider t > 0:
  shape = matrix(c(1,0,1/2,sqrt(3)/2,0,0), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = 0, to = 1, length.out = 1001)
  X = input_operation(shape, seq_grid, "homothety")
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 2, nb_text = 5, asp = 1)
  
  # We can plot all 6 coordinates of the map.
  # It looks simple functions, however it is difficult to get simple
  # equations of them.
  plot(seq_grid, mapX[1,1,], type = "l")
  plot(seq_grid, mapX[2,1,], type = "l")
  plot(seq_grid, mapX[3,1,], type = "l")
  plot(seq_grid, mapX[1,2,], type = "l")
  plot(seq_grid, mapX[2,2,], type = "l")
  plot(seq_grid, mapX[3,2,], type = "l")
  
  # With calculations, we obtain (assuming t>0):
  # x4 = (1, 0)
  #
  # x5 = (\frac{1}{2*\sqrt{a}}) (2*\sqrt{a}-2+t, \sqrt{3}*t)
  # by defining a = t^2 -t + 1
  #
  # x6 = (\frac{sqrt{b}-1}{sqrt{b}}) (b/2, \frac{\sqrt{3}}{2}*\frac{t}{\sqrt{a}}) 
  # by defining b = \frac{1}{alpha}(2*\sqrt{a} - 2 + t)
  #
  # This becomes tedious! And even more after. Not the correct point of view.
  #
  # With t --> 0, we have: x5 = (0, 0), x6 = (0, -1).
  #
  # With t --> infty, we have: 
  # x5 = (1/2)(3, \sqrt{3}), 
  # x6 = (1/2)(3-\sqrt{3}, \sqrt{3}-1)
  # 
  # Some plots to represent x5 and x6:
  t = seq(from = 0, to = 10, length.out = 100)
  a = t^2 - t + 1
  b = (1/sqrt(a))*(2*sqrt(a)-2+t)
  
  x5 = cbind(b/2, 
             (sqrt(3)/2)*(t/sqrt(a)))
  plot(x5, xlab = "dim 1", ylab = "dim 2", main = "x5(t)")
  text(x5[,1], x5[,2], round(t, 1), cex = 0.7)
  
  x6 = cbind((1/2)*(b/sqrt(b)*(sqrt(b)-1)), 
             (sqrt(3)/2)*(b/sqrt(b)*(sqrt(b)-1))/sqrt(b) )
  plot(x6, xlab = "dim 1", ylab = "dim 2", main = "x6(t)")
  text(x6[,1], x6[,2], round(t, 1), cex = 0.7)
}

######################
# Plots examples (3) # (input_random)
######################
# Examples of plots with input_random and mapX_func
if(verbose) {
  ##
  # Example 1 (triangle, 2 dimensional space, 3 moving value) 
  ##
  # All NA here. Difficult to see anything on the plot with all the symmetries
  shape = matrix(c(NA,NA,NA,NA,NA,NA), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 10000, 0, 1)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p")
  
  ##
  # Example 2 (triangle, 2 dimensional space, 3 moving value) 
  ##
  # Five NA. One less symmetry, but still difficult to see anything.
  shape = matrix(c(NA,NA,NA,NA,NA,0), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 10000, 0, 1)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p")
  
  ##
  # Example 3 (triangle, 2 dimensional space, 3 moving value) 
  ##
  # Four NA. One less symmetry, looks strange.
  shape = matrix(c(NA,NA,NA,NA,0,0), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 10000, 0, 1)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p")
  
  ##
  # Example 4 (triangle, 2 dimensional space, 3 moving value) !!
  ##
  # Three NA. Looks great!
  shape = matrix(c(NA,NA,NA,0,0,0), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 10000, 0, 1)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p")
  
  ##
  # Example 4 (triangle, 2 dimensional space, 3 moving value)
  ##
  # Two NA. Here we observe that it only remains one degree of freedom.
  shape = matrix(c(NA,NA,1,0,0,0), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 10000, 0, 1)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p")
  
  ##
  # Example 5 (triangle, 2 dimensional space, 3 moving value)
  ##
  # One NA. Still one degree of freedom.
  shape = matrix(c(NA,3,1,0,0,0), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 1000, 0, 1)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p")

  ##
  # Example 6 (triangle, 2 dimensional space, 3 moving value) 
  ##
  # There are some symmetries
  shape = matrix(c(0,0,NA,0,NA,NA), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 10000, -1, 1)
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p")
  
  ##
  # Example 7 (triangle, 2 dimensional space, 3 moving value) 
  ##
  # Like example 6, but removing one symmetry. Quite beautiful.
  shape = matrix(c(0,0,NA,0,NA,NA), nrow = 3, ncol = 2, byrow = TRUE)
  X = input_random(shape, 10000, 0, 1) # !!
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, asp = 1, type = "p")
}

#############
# Plots (4) # (max norm and one norm)
#############
# Define maps for the infinite norm and one norm
map_max = map_custom(type_norm = max_norm,
                     epsilon = 1e-5,
                     nwalk = 1000)

map_one = map_custom(type_norm = one_norm,
                     epsilon = 1e-5,
                     nwalk = 1000)

# Norm and rotations
if(verbose) {
  # Contrary to the Euclidian norm, the mapping does not verify R*map = map*R
  # Since the unit ball looks like a square, we obtain this kind of symmetry
  # in the output.
  
  ##
  # Max norm and rotation
  ##
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE) # !!
  seq_grid = seq(from = -2, to = 5, length.out = 5001)
  X = input_operation(shape, seq_grid, "rotation")
  mapX = mapX_func(X, map_func = map_max)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 2, nb_text = 5, asp = 1)
  
  ##
  # One norm and rotation
  ##
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -2, to = 5, length.out = 5001)
  X = input_operation(shape, seq_grid, "rotation")
  mapX = mapX_func(X, map_func = map_one)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 2, nb_text = 5, asp = 1)
  
  ##
  # Euclidian norm and rotation
  ##
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = -2, to = 5, length.out = 5001)
  X = input_operation(shape, seq_grid, "rotation")
  mapX = mapX_func(X)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 2, nb_text = 5, asp = 1)
}

# Norm and homotheties
if(verbose) {
  ##
  # Example 1: One norm
  ##
  shape = matrix(c(1,0,1/2,sqrt(3)/2,0,0), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = 0, to = 10, length.out = 2001)
  X = input_operation(shape, seq_grid, "homothety")
  mapX = mapX_func(X, map_func = map_one)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 2, nb_text = 5, asp = 1)
  
  ##
  # Example 2: max norm
  ##
  shape = matrix(c(1,0,1/2,sqrt(3)/2,0,0), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = 0, to = 100, length.out = 2001)
  X = input_operation(shape, seq_grid, "homothety")
  mapX = mapX_func(X, map_func = map_max)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 2, nb_text = 5, asp = 1)
  
  ##
  # Example 2: max norm, other shape
  ##
  shape = matrix(c(0,0,1,0,1,1), nrow = 3, ncol = 2, byrow = TRUE)
  seq_grid = seq(from = 0, to = 3, length.out = 501)
  X = input_operation(shape, seq_grid, "homothety")
  mapX = mapX_func(X, map_func = map_max)
  plot_mapX(seq_grid, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
            nb_polygon = 2, nb_text = 5, asp = 1)
  
  # Some coordinates seems easy to represent, others are cut into many pieces.
  plot(seq_grid, mapX[1,1,], type = "l") # ??
  plot(seq_grid, mapX[2,1,], type = "l") # ??
  plot(seq_grid, mapX[3,1,], type = "l") # ??
  plot(seq_grid, mapX[1,2,], type = "l") # neat
  plot(seq_grid, mapX[2,2,], type = "l") # neat
  plot(seq_grid, mapX[3,2,], type = "l") # neat
}

#########
# Study # (len = 3, m = 2, rule_shift, norm = Euclidian)
#########
# Mathematically, we can restrict the study from 6 to 1 dimension when
# we consider rule_shift, the Euclidian norm (default values),
# and a shape being a triangle in dimension 2
# 
# We then need to study those initial points, for t in ]0, pi]
# (t < 0 is symmetric).
# x1 = (0, 0)
# x2 = (1, 0)
# x3 = (cos(t), sin(t)).
#
# For this purpose, we define the following array X:
shape = matrix(c(0,0,1,0,NA,NA), nrow = 3, ncol = 2, byrow = TRUE)
t = seq(from = 0.0001, to = pi, length.out = 5001) # this is seq_grid
X = array(rep(shape, length(t)), dim = c(3, 2, length(t)))
z = exp(1i*t)
X[3,1,] = Re(z)
X[3,2,] = Im(z)
# plot(X[3,1,], X[3,2,], asp = 1, type = "l") # this is a half-circle

# This is the map we want to study
mapX = mapX_func(X)
plot_mapX(t, mapX, xaxis = 1, yaxis = 2, cex = 0.1,
          nb_polygon = 7, nb_text = 7, asp = 1)

# Some numerical checks
if(verbose) {
  ##
  # Mean value
  ##
  # Mean value of the triangle outputs. It looks like a arc of circle,
  # but it should not be this.
  mean_map = apply(mapX[,,], c(2,3), mean)
  plot(mean_map[1,], mean_map[2,], asp = 1, type = "l")
  
  ##
  # Plot of each output coordinates.
  ##
  # The functions seems simple...
  plot(t, mapX[1,1,], type = "l")
  plot(t, mapX[2,1,], type = "l")
  plot(t, mapX[3,1,], type = "l")
  plot(t, mapX[1,2,], type = "l")
  plot(t, mapX[2,2,], type = "l")
  plot(t, mapX[3,2,], type = "l")
  
  # For example, mx_2 looks like a circle
  x = mapX[2,1,]
  y = mapX[2,2,]
  plot(x,y)
  
  # We can try to estimate the center of the 'circle'
  grid = expand.grid(seq(from = 0.95, to = 1.05, length.out = 100),
                     seq(from = -0.30, to = -0.20, length.out = 100))
  grid = cbind(grid, rep(NA, nrow(grid)))
  for(i in 1:nrow(grid)) {
    a = grid[i,1]
    b = grid[i,2]
    estim_r = sqrt((x-a)^2 + (y-b)^2)
    grid[i,3] = sd(estim_r)
  }
  i = which.min(grid[,3])
  center_x = grid[i,1] # ~1.003
  center_y = grid[i,2] # ~0.256
  estim_r = sqrt((x-center_x)^2 + (y-center_y)^2)
  
  # The approximation is good
  r = mean(estim_r) # ~0.255
  center = c(a, b)
  
  # However it's not really a circle. For example, we can see some variations
  # in the radius:
  plot(t, estim_r, xlab = "t", 
       ylab = "Distance from 'center' to element t")
  
  ##
  # All 2D cross-plots
  ##
  plot(mapX[1,1,], mapX[2,1,], type = "l", asp = 1)
  plot(mapX[1,1,], mapX[3,1,], type = "l", asp = 1)
  plot(mapX[1,1,], mapX[1,2,], type = "l", asp = 1)
  plot(mapX[1,1,], mapX[2,2,], type = "l", asp = 1)
  plot(mapX[1,1,], mapX[3,2,], type = "l", asp = 1)
  #
  plot(mapX[2,1,], mapX[3,1,], type = "l", asp = 1)
  plot(mapX[2,1,], mapX[1,2,], type = "l", asp = 1)
  plot(mapX[2,1,], mapX[2,2,], type = "l", asp = 1)
  plot(mapX[2,1,], mapX[3,2,], type = "l", asp = 1)
  #
  plot(mapX[3,1,], mapX[1,2,], type = "l", asp = 1)
  plot(mapX[3,1,], mapX[2,2,], type = "l", asp = 1)
  plot(mapX[3,1,], mapX[3,2,], type = "l", asp = 1)
  #
  plot(mapX[1,2,], mapX[2,2,], type = "l", asp = 1)
  plot(mapX[1,2,], mapX[3,2,], type = "l", asp = 1)
  #
  plot(mapX[2,2,], mapX[3,2,], type = "l", asp = 1)
}
