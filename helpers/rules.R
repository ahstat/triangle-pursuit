##
# Define the rule of the sequence to compute x[i] from x[1:(i-1)]
##

# Euclidian distance
euc_norm = function(vect) {
  return(norm(vect, type="2"))
}

# Infinite distance
max_norm = function(vect) {
  return(max(abs(vect)))
}

# 1-norm distance
one_norm = function(vect) {
  return(sum(abs(vect)))
}

# modulo function with positive result
`%mod%` = function(i, len) {
  modulo = i %% len
  modulo = ifelse(modulo == 0, modulo+len, modulo)
  return(modulo)
}

# Walk a distance 1 from x to y
# (return NaN if x==y)
go_from_to = function(x, y, type_norm = euc_norm) {
  vect = x-y
  out = x - vect / type_norm(vect)
  return(out)
}

# Rule shift:
# x[i] defined as the step from x[i-1] to the direction x[i-len] of length 1
#
#  Example: x1, x2, x3 initial points. 
#    Then x4 from x3 to x1 with length 1,
#         x5 from x4 to x2 with length 1,
#         x6 from x5 to x3 with length 1,
#         x7 from x6 to x4 with length 1,
#         xn from x(n-1) to x(n-3) with length 1.
rule_shift = function(x, i, len, type_norm) {
  out = go_from_to(x[i-1,], x[i-len,], type_norm)
  return(out)
}

# Rule origin:
# x[i] defined as the step from from x[i-1] to the direction x[modulo] of
#   length 1. Here, modulo is one of the first elements of sequence x,
#   from x[1] to x[len], and such that: modulo = i (mod len).
#
#  Example: x1, x2, x3 initial points. 
#    Then x4 from x3 to x1 with length 1,
#         x5 from x4 to x2 with length 1,
#         x6 from x5 to x3 with length 1,
#         x7 from x6 to x1 with length 1,
#         xn from x(n-1) to x(n mod 3) with length 1.
rule_origin = function(x, i, len, type_norm) {
  modulo = i %mod% len
  out = go_from_to(x[i-1,], x[modulo,], type_norm)
  return(out)
}