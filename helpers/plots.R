# Take the dimension xaxis is xaxis (from 1 to dim(mapX)[2])
# NA for grid as axis
get_feat = function(mapX, seq_grid, axis) {
  if(is.na(axis)) { # take the grid as x_features
    nb_points = dim(mapX)[1]
    feats = matrix(rep(seq_grid, nb_points), nrow = nb_points, byrow = TRUE)
    rownames(feats) = rownames(mapX)
    colnames(feats) = dimnames(mapX)[[3]]
  } else { # of the form "dim i"
    feats = mapX[,axis,]
  }
  return(feats)
}

get_lab = function(axis) {
  if(is.na(axis)) {
    return("grid values")
  } else {
    return(paste("dim", axis))
  }
}

plot2 = function(..., add = TRUE) {
  if(add == FALSE) {
    plot(...)
  } else {
    lines(...)
  }
}

# xaxis and yaxis should be in 1:dim(mapX)[2], or to be NA.
# If NA, then the grid will be taken,
# Else, the corresponding dimension will be taken
plot_mapX = function(seq_grid, mapX, 
                     xaxis = NA, yaxis = 1,
                     cols = NA, nb_text = NA, nb_polygon = NA,
                     main = NA, cex = 0.8, type = "p", 
                     fade = TRUE, ...) {
  if(is.na(cols)[1]) {
    cols = c("red", "blue", "darkgreen", "orange", "purple1",
             "black", "yellow2", "turquoise2", "springgreen1", "snow3")
  }
  
  if(length(dim(mapX)) != 3) {
    stop("Plot with only one moving variable.")
  }
  
  # x feats
  x_feats = get_feat(mapX, seq_grid, xaxis)
  xlim = range(x_feats, na.rm=TRUE)
  xlab = get_lab(xaxis)
  
  # y feats
  y_feats = get_feat(mapX, seq_grid, yaxis)
  ylim = range(y_feats, na.rm=TRUE)
  ylab = get_lab(yaxis)
  
  dims = dim(x_feats)[1]
  for(i in 1:dims) {
    colfunc_i = colorRampPalette(c(cols[i],"lightgray"))
    add = ifelse(i == 1, FALSE, TRUE)
    size = ifelse(fade, ncol(x_feats), 1)
    plot2(x_feats[i,], y_feats[i,], type = type, 
          ylim = ylim, xlim = xlim,
          ylab = ylab, xlab = xlab, 
          col = colfunc_i(size), pch = 19, cex = cex,
          main = main,
          add = add, ...)
    
    # Add the grid number on the plot
    if(!is.na(nb_text)) {
      idx = seq(from = 1, to = length(seq_grid), length.out = nb_text)
      text(x_feats[i,][idx], y_feats[i,][idx], round(seq_grid[idx], 2),
           col = colfunc_i(1), pos = 2, cex = 0.7)
    }
  }
  
  if(!is.na(nb_polygon)) {
    colfunc_polygon = colorRampPalette(c("black","lightgray"))
    col_polygon = colfunc_polygon(nb_polygon)
    idxs = seq(from = 1, to = length(seq_grid), length.out = nb_polygon)
    for(idx in 1:length(idxs)) {
      j = idxs[idx]
      polygon(x_feats[,j], y_feats[,j],
              col = col_polygon[idx], density = 0)
    }
  }
}