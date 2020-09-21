# Default par() settings for Voronoi tessellations 
.vis_polytope_setup = function(VOR) {
  
  ## Determine plotting bounds 
  xlim = c(VOR$lb[1], VOR$ub[1])
  ylim = c(VOR$lb[2], VOR$ub[2])
  
  ## Plot a blank window 
  par(xaxs = 'i', yaxs = 'i', mar = c(4.1, 4.1, 1.1, 1.1))
  plot(0, cex=0, xlab = 'x', ylab = 'y', xlim = xlim, ylim = ylim)
  
  ## Store par for later use
  VOR$set_vis_par(par(no.readonly = T))
}


#' Tableau 20 Color Palette 
#' @return a named vector containing Tableau 20 colors 
tableau20 = function(col = NULL) {
  mycolors = c("blue" = "#1F77B4", "lightblue" = "#73C2FB", 
               "orange" = "#ff7f0e", "lightorange" = "#ffd27f", 
               "green" = "#2ca02c", "lightgreen" = "#98fb98", 
               "red" = "#d62728", "lightred" = "#ffb09c", 
               "purple" = "#b660cd", "lightpurple" = "#e4a0f7", 
               "yellow" = "#ffdb58", "lightyellow" = "#fdfd96", 
               "teal" = "#17becf", "lightteal" = "#c8ffff", 
               "gray" = "#88807b", "lightgray" = "#c7c6c1", 
               "brown" = "#8c564b", "lightbrown" = "#ceb180", 
               "pink" = "#ff6fff", "lightpink" = "#fde6fa")
  
  if(is.null(col)) {
    return(mycolors)
  } else {
    return(mycolors[col])
  }
}



#' Compute the convex hull of vertices of a polytope
#'
#' Polytope definition should be given in the H-representation Ax <= b. 
#' This function uses Komei Fukuda's cdd library to convert it to a V-representation. 
#'
#' @param A the LHS constraint matrix 
#' @param b the RHS constraint vector 
#' @return a list with components: 
#' \enumerate{
#' \item V a matrix with vertex points in its rows
#' \item Vpath a 2 column edge list matrix whose rows specify vertex indices (index to rows of V)
#'   that are connected in the convex hull path. Note that this path has its redundancies removed 
#'   (i.e., a connection between vertex 1 & 2 is given only by a row (1,2), not by (2,1))
#'}
#' @export
polytope_chull = function(A, b) {
  
  ## Make the H-polytope representation in rcdd
  Hpoly = rcdd::makeH(a1 = A, b1 = b)
  
  ## Convert it to a V-representation with rcdd
  Vpoly = rcdd::scdd(input = Hpoly, adjacency = T)
  
  ## Now we need to parse the adjacency list returned in Vpoly 
  ## Create a temporary function to do this
  parse_rcdd_adjacency = function(i) {
    tmp = Vpoly$adjacency[[i]]
    tmp=tmp[tmp>i]
    if(length(tmp)==0) return(NULL)
    return(cbind(i, tmp))
  }
  
  ## Apply the function 
  vpath = sapply(1:length(Vpoly$adjacency), parse_rcdd_adjacency)
  vpath = do.call(rbind, vpath)
  colnames(vpath) = NULL 
  
  
  ## Decode the path to get the connecting order 
  Vpathtmp = vpath
  Vorder = Vpathtmp[1,]
  Vpathtmp = Vpathtmp[-1,,drop=F]
  
  while(nrow(Vpathtmp) > 1) {
    ij = which(Vpathtmp == Vorder[1], arr.ind = T)
    if(ij[2]==1)
      Vorder = c(Vpathtmp[ij[1],2], Vorder)
    else
      Vorder = c(Vpathtmp[ij[1],1], Vorder)
    
    Vpathtmp = Vpathtmp[-ij[1],,drop=F]
  }
  
  
  ## Return our outputs
  return(list(V = Vpoly$output[,3:ncol(Vpoly$output)], Vpath = vpath, Vorder = Vorder))
}



#' Visualize polytope constraints
#'
#' The boundaries x : Ax = b will be plotted as lines
#'
#' @param A the constraint matrix (must be 2-dimensional)
#' @param b the constraint RHS
#' @param add boolean, whether to generate a new plot or add to existing
#' @param edge.lty edge linetype
#' @param edge.lwd edge linewidth
#' @param edge.col edge color
#' @return nothing, a plot will be generated
#' @details The plotting parameters *.lty, *.lwd and *.col are passed to R's plot() command. \cr
#' @export
vis_polytope_constraints = function(A, b, add=F, edge.lty = 1, edge.lwd = 1, edge.col='magenta') {
  
  ## Setup plot, if requested
  if(!add) {
    ## Determine the vertices of the polytope.
    ## This is needed to compute the plotting bounds
    #verts = polytope_vertices_2d(A, b)
    verts = polytope_chull(A = A, b = b)$V
    
    xlim = range(verts[,1])
    ylim = range(verts[,2])
    
    opar = par(no.readonly = T)
    on.exit(par(opar))
    par(xaxs = 'i', yaxs = 'i', mar = c(4.1, 4.1, 1.1, 1.1))
    plot(0, cex=0, xlim=xlim, ylim=ylim, asp=1, xlab = 'x', ylab = 'y')
  }
  
  ## Loop over each constraint
  for(i in 1:nrow(A)) {
    
    ## If both x & y coefficients are nonzero, plot the usual slope-intercept form of a line
    if(A[i,1] != 0 && A[i,2] != 0) {
      intercept = b[i] / A[i,2]
      slope = -A[i,1] / A[i,2]
      abline(a = intercept, b = slope, lty=edge.lty, lwd=edge.lwd, col=edge.col)
      
      ## If only the x coefficient is zero, plot a horizontal line
    } else if(A[i,1] == 0) {
      abline(h = b[i] / A[i,2], lty=edge.lty, lwd=edge.lwd, col=edge.col)
      
      ## If only the y coefficient is zero, plot a vertical line
    } else if(A[i,2] == 0) {
      abline(v = b[i] / A[i,1], lty=edge.lty, lwd=edge.lwd, col=edge.col)
    }
  } ## end loop over constraints
}



#' Visualize a (bounded) polytope
#'
#' The boundary region of a bounded polytope will be plotted as a set of intersecting line segments.
#' @param A the constraint matrix (must be 2-dimensional)
#' @param b the constraint RHS
#' @param center a flag denoting the type of center where the polytope label will be plotted as text. Can be one of
#' \enumerate{
#' \item NULL (default), no center defined
#' \item a vector (length = 2) giving the (x,y) coords of the center
#' \item "chebyshev", in which case the Chebyshev center will be calculated and used 
#' \item "analytic", in which case the analytic center will be calculated and used
#' }
#' @param add boolean, whether to generate a new plot or add to existing
#' @param fill.col color to fill polygon. Default = NULL means no fill used. 
#' @param edge.lty edge linetype
#' @param edge.lwd edge linewidth
#' @param edge.col edge color
#' @param vertex.pch the vertex point style
#' @param vertex.cex the vertex point size, set = 0 to suppress vertex point plotting
#' @param vertex.col the vertex point color
#' @param center.pch the center point style
#' @param center.cex the center point size, set = 0 to suppress center point plotting if \code{center} != NULL
#' @param center.col the center point color
#' @param text a string (or numeric value, which will be cast as string) giving the label text which is plotted at (x,y) coords defined by \code{center}.
#' Default = NULL (nothing).
#' @param text.cex the text label size
#' @param text.col the text label color
#' @return nothing, a plot is generated
#' @details \code{center} must be given as something OTHER than NULL if either a point or label is requested to be plotted at the polytope center.
#' If only a text label at the center is desired, set center.cex = 0.
#' If only a point at the center is desired, set label = NULL
#' @export
vis_polytope = function(A, b, center = NULL, add = F,
                        fill.col = NULL, 
                        edge.lty = 1, edge.lwd = 1, edge.col='navy',
                        vertex.pch = 16, vertex.cex = 0.75, vertex.col = 'navy',
                        center.pch = 16, center.cex = 0.75, center.col = 'navy',
                        text = NULL, text.cex = 0.75, text.col = 'navy') {
  
  ## *** Checks:
  ## Make sure A is just 2-d
  if(ncol(A) != 2) stop("Polytope must be 2-dimensional.")
  ## Make sure constraint system dimensions agree
  if(nrow(A) != length(b)) stop("nrow(A) must = length(b).")
  
  ## Process the type of center point requested 
  ## If center is given, it must be either a length=2 vector defining the (x,y) coords of the center, 
  ## or the string "analytic" in which case the analytic center will be calculated and plotted 
  if(is.numeric(center)) {
    ## Check dimensions 
    center = c(center)
    if(length(center)!=2) stop("If given as numeric, center must be a vector of length=2 giving (x,y) coords of polytope center.")
  } else if(isTRUE(tolower(center) == "analytic")) {
    ## Compute analytic center 
    center = c(VorVQ::polytope_analytic_center(A = A, b = b))
  } else if(isTRUE(tolower(center) == "chebyshev")) {
    center = c(VorVQ::polytope_chebyshev_center(A = A, b = b))
  } else if(is.null(center)) {
    ## Do nothing, already marked flags as F
  } else {
    stop("If center is given, it must be a vector of length=2, or one of 'chebyshev' or 'analytic'.")
  }
  
 
  ## Get the convex hull of the polytope
  ## This returns a list with components "V" and "Vpath" giving the vertices and edges connecting them in the convex hull 
  poly_hull = VorVQ::polytope_chull(A = A, b = b)
  
  
  ## Setup plot, if requested
  if(!add) {
    ## Compute the range of the vertices, +/- 2% 
    xlim = range(poly_hull$V[,1])
    xlim = xlim + c(-.02*diff(xlim), .02*diff(xlim)) ## add a buffer range
    ylim = range(poly_hull$V[,2])
    ylim = ylim + c(-.02*diff(ylim), .02*diff(ylim)) ## add a buffer range
    
    opar = par(no.readonly = T)
    on.exit(par(opar))
    par(xaxs = 'i', yaxs = 'i', mar = c(4.1, 4.1, 1.1, 1.1))
    plot(0, cex=0, xlim=xlim, ylim=ylim, asp=1, xlab = 'x', ylab = 'y')
  }
  
 
  ## Fill the polygon, if requested 
  if(!is.null(fill.col)) {
    polygon(poly_hull$V[poly_hull$Vorder,,drop=F], border = NA, col = fill.col)
  }
  
  ## Draw the edges of the polytope
  if(edge.lwd > 0) {
    segments(x0 = poly_hull$V[poly_hull$Vpath[,1],1],
             y0 = poly_hull$V[poly_hull$Vpath[,1],2],
             x1 = poly_hull$V[poly_hull$Vpath[,2],1],
             y1 = poly_hull$V[poly_hull$Vpath[,2],2],
             lty = edge.lty, lwd = edge.lwd, col = edge.col)    
  }
  
  
  ## Add the vertices, if requested
  if(vertex.cex > 0) {
    points(x = poly_hull$V[,1], y = poly_hull$V[,2], pch = vertex.pch, cex = vertex.cex, col = vertex.col)
  }
  
  ## Add the center point, if requested
  if(!is.null(center) && center.cex > 0) {
    points(x = center[1], y = center[2], pch = center.pch, cex = center.cex, col = center.col)
  }
  
  ## Add the labels, if requested
  if(!is.null(center) && !is.null(text) && label.cex > 0) {
    text = as.character(text)
    text(x = center[1], y = center[2], labels = text, cex = text.cex, col = text.col)
  }
  
}



#' Visualize the first-order Voronoi cells generated by prototypes
#'
#' @param VOR a Voronoi object, instantiated by \code{VOR$new()}
#' @param add boolean, whether to generate a new plot or add to existing
#' @param edge.lty edge linetype
#' @param edge.lwd edge linewidth
#' @param edge.col edge color
#' @param vertex.pch the vertex point style
#' @param vertex.cex the vertex point size, set = 0 to suppress vertex point plotting
#' @param vertex.col the vertex point color
#' @return nothing, a plot is generated
#' @export 
vis_vor1_polytopes = function(VOR, add = F, 
                              edge.lty = 1, edge.lwd = 1, edge.col=tableau20("gray"),
                              vertex.pch = 16, vertex.cex = 0.75, vertex.col = tableau20("gray")) {
  
  ## *** Checks:
  # Make sure the input is the proper class 
  if(!class(VOR)=="Rcpp_VOR") stop("Input VOR must be an object of class VOR")
  # Make sure the DADJ has been calculated. If not, do so 
  #if(nrow(VOR$DADJ)==0) stop("Must call $calc_DADJ before any visualization")
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  ## *** Setup plot, if requested
  if(add) {
    par(VOR$vis_par)
  } else {
    .vis_polytope_setup(VOR)
  }
 
  
  ## Build a new set of plotting bounds to extend the polytope vertices 
  # at the border past the plot window 
  newlb = rep(min(VOR$lb),2)
  newub = rep(max(VOR$ub),2)
  newrng = newub - newlb
  newlb = newlb - 0.05*newrng
  newub = newub + 0.05*newrng 
  
  ## *** Plot, loop over each cell 
  for(i in VOR$vor1_active) {
    if(nrow(VOR$DADJ)>0) {
      poly = VorVQ::vor1_polytope(W = VOR$W, iidx = VOR$vor1_active[i], DADJ = VOR$DADJ, lb = newlb, ub = newub, rmv_redundancies = F)  
    } else {
      poly = VorVQ::vor1_polytope(W = VOR$W, iidx = VOR$vor1_active[i], lb = newlb, ub = newub, rmv_redundancies = F)
    }
    
    VorVQ::vis_polytope(A = poly$A, b = poly$b, center = NULL, add = T, 
                           fill.col = NULL, 
                           edge.lty = edge.lty, edge.lwd = edge.lwd, edge.col = edge.col, 
                           vertex.pch = vertex.pch, vertex.cex = vertex.cex, vertex.col = vertex.col)
  }
  
 
}


#' Visualize the second-order Voronoi cells generated by prototypes
#'
#' @param VOR a Voronoi object, instantiated by \code{VOR$new()}
#' @param add boolean, whether to generate a new plot or add to existing
#' @param edge.lty edge linetype
#' @param edge.lwd edge linewidth
#' @param edge.col edge color
#' @param vertex.pch the vertex point style
#' @param vertex.cex the vertex point size, set = 0 to suppress vertex point plotting
#' @param vertex.col the vertex point color
#' @return nothing, a plot is generated
#' @export 
vis_vor2_polytopes = function(VOR, add = F, 
                              edge.lty = 1, edge.lwd = 1, edge.col=tableau20("lightgray"),
                              vertex.pch = 16, vertex.cex = 0.75, vertex.col = tableau20("lightgray")) {
  
  ## *** Checks:
  # Make sure the input is the proper class 
  if(!class(VOR)=="Rcpp_VOR") stop("Input VOR must be an object of class VOR")
  # Make sure the DADJ has been calculated. If not, do so 
  #if(nrow(VOR$DADJ)==0) stop("Must call $calc_DADJ before any visualization")
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  ## *** Setup plot, if requested
  ## *** Setup plot, if requested
  if(add) {
    par(VOR$vis_par)
  } else {
    .vis_polytope_setup(VOR)
  }
  
  
  ## Build a new set of plotting bounds to extend the polytope vertices 
  # at the border past the plot window 
  newlb = rep(min(VOR$lb),2)
  newub = rep(max(VOR$ub),2)
  newrng = newub - newlb
  newlb = newlb - 0.05*newrng
  newub = newub + 0.05*newrng 
  
  
  ## *** Plot, loop over each cell 
  for(i in 1:nrow(VOR$vor2_active)) {
    if(nrow(VOR$DADJ) > 0) {
      poly = VorVQ::vor2_polytope(W = VOR$W, iidx = VOR$vor2_active[i,1], jidx = VOR$vor2_active[i,2], DADJ = VOR$DADJ, lb =newlb, ub = newub, rmv_redundancies = F)  
    } else {
      poly = VorVQ::vor2_polytope(W = VOR$W, iidx = VOR$vor2_active[i,1], jidx = VOR$vor2_active[i,2], lb = newlb, ub = newub, rmv_redundancies = F)  
    }
    
    VorVQ::vis_polytope(A = poly$A, b = poly$b, center = NULL, add = T, 
                           fill.col = NULL, 
                           edge.lty = edge.lty, edge.lwd = edge.lwd, edge.col = edge.col, 
                           vertex.pch = vertex.pch, vertex.cex = vertex.cex, vertex.col = vertex.col)
  }

}



#' Visualize the first-order Voronoi MVIEs
#'
#' @param VOR a Voronoi object, instantiated by \code{VOR$new()}
#' @param add boolean, whether to generate a new plot or add to existing
#' @param col color for plotting ellipse and its center
#' @param ell.lwd linewidth of ellipse 
#' @param center.pch the center point style
#' @param center.cex the center point size, set = 0 to suppress center point plotting
#' @return nothing, a plot is generated
#' @export 
vis_vor1_MVIE = function(VOR, add = F,
                         col = tableau20("blue"), 
                         ell.lwd = 1.5, center.pch = 16, center.cex = 0.75) {
  
  ## *** Checks:
  # Make sure the input is the proper class 
  if(!class(VOR)=="Rcpp_VOR") stop("Input VOR must be an object of class VOR")
  # Make sure the DADJ has been calculated. If not, do so 
  if(nrow(VOR$vor1_MVIE_c)!=length(VOR$vor1_active)) stop("Must call $calc_vor1_MVIE before any visualization")
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  ## *** Setup plot, if requested
  if(add) {
    par(VOR$vis_par)
  } else {
    .vis_polytope_setup(VOR)
  }
  
  
  ## *** Plot, loop over each cell 
  for(i in 1:length(VOR$vor1_active)) {
    car::ellipse(center = c(VOR$vor1_MVIE_c[i,]), shape = VOR$vor1_MVIE_E[,,i], radius = 1, draw = T, col = col, lwd = ell.lwd, center.pch = center.pch, center.cex = center.cex, add = T)
  }

}


#' Visualize the second-order Voronoi MVIEs
#'
#' @param VOR a Voronoi object, instantiated by \code{VOR$new()}
#' @param add boolean, whether to generate a new plot or add to existing
#' @param col color for plotting ellipse and its center
#' @param ell.lwd linewidth of ellipse 
#' @param center.pch the center point style
#' @param center.cex the center point size, set = 0 to suppress center point plotting
#' @return nothing, a plot is generated
#' @export 
vis_vor2_MVIE = function(VOR, add = F,
                         col = tableau20("lightblue"), 
                         ell.lwd = 1.5, center.pch = 16, center.cex = 0.75) {
  
  ## *** Checks:
  # Make sure the input is the proper class 
  if(!class(VOR)=="Rcpp_VOR") stop("Input VOR must be an object of class VOR")
  # Make sure the DADJ has been calculated. If not, do so 
  if(nrow(VOR$vor2_MVIE_c)!=nrow(VOR$vor2_active)) stop("Must call $calc_vor2_MVIE before any visualization")
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  ## *** Setup plot, if requested
  if(add) {
    par(VOR$vis_par)
  } else {
    .vis_polytope_setup(VOR)
  }
  
  
  ## *** Plot, loop over each cell 
  for(i in 1:nrow(VOR$vor2_active)) {
    car::ellipse(center = c(VOR$vor2_MVIE_c[i,]), shape = VOR$vor2_MVIE_E[,,i], radius = 1, draw = T, col = col, lwd = ell.lwd, center.pch = center.pch, center.cex = center.cex, add = T)
  }
  
}




#' Visualize the first-order Voronoi Dikin Ellipses
#'
#' @param VOR a Voronoi object, instantiated by \code{VOR$new()}
#' @param add boolean, whether to generate a new plot or add to existing
#' @param col color for plotting ellipse and its center
#' @param ell.lwd linewidth of ellipse 
#' @param center.pch the center point style
#' @param center.cex the center point size, set = 0 to suppress center point plotting
#' @return nothing, a plot is generated
#' @export 
vis_vor1_Dikin = function(VOR, add = F,
                         col = tableau20("green"), 
                         ell.lwd = 1.5, center.pch = 16, center.cex = 0.75) {
  
  ## *** Checks:
  # Make sure the input is the proper class 
  if(!class(VOR)=="Rcpp_VOR") stop("Input VOR must be an object of class VOR")
  # Make sure the Dikin ellipse has been calculated. 
  if(dim(VOR$vor1_Dikin_E)[3]!=length(VOR$vor1_active)) stop("Must call $calc_vor1_Dikin before any visualization")
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  ## *** Setup plot, if requested
  if(add) {
    par(VOR$vis_par)
  } else {
    .vis_polytope_setup(VOR)
  }
  
  
  ## *** Plot, loop over each cell 
  centers = vor$get_vor1_centers()
  for(i in 1:length(VOR$vor1_active)) {
    car::ellipse(center = centers[i,], shape = VOR$vor1_Dikin_E[,,i], radius = 1, draw = T, col = col, lwd = ell.lwd, center.pch = center.pch, center.cex = center.cex, add = T)
  }
  
}



#' Visualize the second-order Voronoi Dikin Ellipses
#'
#' @param VOR a Voronoi object, instantiated by \code{VOR$new()}
#' @param add boolean, whether to generate a new plot or add to existing
#' @param col color for plotting ellipse and its center
#' @param ell.lwd linewidth of ellipse 
#' @param center.pch the center point style
#' @param center.cex the center point size, set = 0 to suppress center point plotting
#' @return nothing, a plot is generated
#' @export 
vis_vor2_Dikin = function(VOR, add = F,
                         col = tableau20("lightgreen"), 
                         ell.lwd = 1.5, center.pch = 16, center.cex = 0.75) {
  
  ## *** Checks:
  # Make sure the input is the proper class 
  if(!class(VOR)=="Rcpp_VOR") stop("Input VOR must be an object of class VOR")
  # Make sure the DADJ has been calculated. If not, do so 
  if(dim(VOR$vor2_Dikin_E)[3]!=nrow(VOR$vor2_active)) stop("Must call $calc_vor2_Dikin before any visualization")
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  ## *** Setup plot, if requested
  if(add) {
    par(VOR$vis_par)
  } else {
    .vis_polytope_setup(VOR)
  }
  
  
  ## *** Plot, loop over each cell
  centers = VOR$get_vor2_centers()
  for(i in 1:nrow(VOR$vor2_active)) {
    car::ellipse(center = centers[i,], shape = VOR$vor2_Dikin_E[,,i], radius = 1, draw = T, col = col, lwd = ell.lwd, center.pch = center.pch, center.cex = center.cex, add = T)
  }
  
}



#' Annotate the first-order Voronoi cells
#'
#' @param VOR a Voronoi object, instantiated by \code{VOR$new()}
#' @param add boolean, whether to generate a new plot or add to existing
#' @param text the text label to be printed in each Voronoi cell. 
#' The default is a single string = 'vorid', which prints the ids in \code{$vor1_active} at the centers of each active first-order cell. 
#' Can also be given as a vector (length = length(\code{$vor1_active})) containing a label for each active first-order cell. 
#' @param text.cex size of plotted text labels 
#' @param text.col color for plotted text labels
#' @param text.font font weight for plotted text labels. Set = 2 for bold. 
#' @details 
#' The text will be plotted at the \code{(x,y)} coordinates given by \code{$get_vor1_centers}. 
#' @return nothing, a plot is generated
#' @export 
vis_vor1_annotate = function(VOR, add = F, 
                             text = "vorid", text.cex = 0.75, text.col = "black", text.font = 1) {
  
  ## *** Checks:
  # Make sure the input is the proper class 
  if(!class(VOR)=="Rcpp_VOR") stop("Input VOR must be an object of class VOR")
  # Get the centers, used to (x,y) coords for plotted labels
  centers = VOR$get_vor1_centers()

  
  ## *** Decode the labels to be printed 
  # If = "vorid", will print a string with "vorid" in each cell 
  # Otherwise, the 'text' vector must have length = length(vor1_active)
  if(any(tolower(na.omit(text)) == "vorid")) {
    text = as.character(VOR$vor1_active)
  } else if(length(text) == length(VOR$vor1_active)) {
    text = as.character(text)
  } else {
    stop("Unknown 'text' input; must be either a single string = 'vorid', or a character vector containing a text label for each cell listed in $vor1_active")
  }
    
  
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  ## *** Setup plot, if requested
  if(add) {
    par(VOR$vis_par)
  } else {
    .vis_polytope_setup(VOR)
  }
  
  
  
  ## Add labels to plot 
  text(centers, labels = text, cex = text.cex, col = text.col, font = text.font)

  
}



#' Annotate the second-order Voronoi cells
#'
#' @param VOR a Voronoi object, instantiated by \code{VOR$new()}
#' @param add boolean, whether to generate a new plot or add to existing
#' @param text the text label to be printed in each Voronoi cell. 
#' The default is a single string = 'vorid', which prints the ids in \code{$vor2_active} at the centers of each active second-order cell. 
#' Can also be given as a vector (length = nrow(\code{$vor2_active})) containing a label for each active second-order cell. 
#' @param text.cex size of plotted text labels 
#' @param text.col color for plotted text labels
#' @param text.font font weight for plotted text labels. Set = 2 for bold. 
#' @details 
#' The text will be plotted at the \code{(x,y)} coordinates given by \code{$get_vor2_centers}. 
#' @return nothing, a plot is generated
#' @export 
vis_vor2_annotate = function(VOR, add = F, 
                             text = "vorid", text.cex = 0.55, text.col = "black", text.font = 1) {
  
  ## *** Checks:
  # Make sure the input is the proper class 
  if(!class(VOR)=="Rcpp_VOR") stop("Input VOR must be an object of class VOR")
  # Try to get the centers, used to (x,y) coords for plotted labels
  # If no centers have been set, this will error 
  centers = VOR$get_vor2_centers() 
  
  
  ## *** Decode the labels to be printed 
  # If = "vorid", will print a string with "vorid" in each cell 
  # Otherwise, the 'text' vector must have length = length(vor1_active)
  if(any(tolower(na.omit(text)) == "vorid")) {
    text = apply(VOR$vor2_active, 1, function(z) paste(z,collapse="-"))
  } else if(length(text) == nrow(VOR$vor2_active)) {
    text = as.character(text)
  } else {
    stop("Unknown 'text' input; must be either a single string = 'vorid', or a character vector containing a text label for each cell listed in $vor2_active")
  }
  
  
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  ## *** Setup plot, if requested
  if(add) {
    par(VOR$vis_par)
  } else {
    .vis_polytope_setup(VOR)
  }
  
  
  ## Add labels to plot 
  text(centers, labels = text, cex = text.cex, col = text.col, font = text.font)
  
  
}
