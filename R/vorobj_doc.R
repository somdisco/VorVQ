#' @name VOR
#' @title The VOR object class
#' 
#' @field W matrix of the prototypes (in rows) which generate the Voronoi tessellation 
#' @field nW the number of prototypes in the tessellation 
#' @field d the dimension of the prototypes \code{= ncols(W)}
#' @field lb a vector of dimension-wise lower bounds for the tessellated region, \code{length = d}
#' @field ub a vector of dimension-wise upper bounds for the tessellated region, \code{length = d}
#' @field CADJ the CADJ adjacency of the learned prototypes
#' @field vor1_active a vector of indices (rows of W) identifying the first-order Voronoi for which approximating ellipsoids will be computed. 
#' Defaults to the active first-order cells in the tessellation (those which contain at least one data observation after the ANN recall, as indicated by CADJ). 
#' @field vor2_active a matrix whose rows contain the (i,j) indices identifying the second-order Voronoi cells for which approximating ellipsoids will be computed. 
#' Defaults to the active second-order cells in the tessellation (those which contain at least one data observation after the ANN recall, as indicated by CADJ). 
#' @field GADJ the Gabriel graph adjacency matrix
#' @field DADJ the Delaunay graph adjacency matrix 
#' @field vor1_centers matrix whose rows store the Chebyshev centers of the first-order Voronoi cells listed in \code{vor1_active}
#' @field vor2_centers matrix whose rows store the Chebyshev centers of the second-order Voronoi cells listed in \code{vor2_active}
#' @field vor1_MVIE_c matrix whose rows store the MVIE centers of the first-order Voronoi cells listed in \code{vor1_active}
#' @field vor1_MVIE_E cube whose slices store the MVIE ellipsoidal rotation matrices of the first-order Voronoi cells listed in \code{vor1_active}
#' @field vor1_MVIE_logdet vector whose elements store the log-determinant of MVIE ellipsoidal rotation matrices of the first-order Voronoi cells listed in \code{vor1_active}
#' @field vor1_MVIE_volratio vector whose elements store the proportional volume of the MVIEs of the first-order Voronoi cells listed in \code{vor1_active}
#' @field vor1_MVIE_status return flag from running the MVIE routine for each of the first-order Voronoi cells listed in \code{vor1_active}
#' @field vor2_MVIE_c matrix whose rows store the MVIE centers of the second-order Voronoi cells listed in \code{vor2_active}
#' @field vor2_MVIE_E cube whose slices store the MVIE ellipsoidal rotation matrices of the second-order Voronoi cells listed in \code{vor2_active}
#' @field vor2_MVIE_logdet vector whose elements store the log-determinant of MVIE ellipsoidal rotation matrices of the second-order Voronoi cells listed in \code{vor2_active}
#' @field vor2_MVIE_volratio vector whose elements store the proportional volume of the MVIEs of the second-order Voronoi cells listed in \code{vor2_active}
#' @field vor2_MVIE_status return flag from running the MVIE routine for each of the second-order Voronoi cells listed in \code{vor2_active}
#' @field vor1_Dikin_E cube whose slices store the Dikin ellipsoidal rotation matrices of the first-order Voronoi cells listed in \code{vor1_active}
#' @field vor2_Dikin_E cube whose slices store the Dikin ellipsoidal rotation matrices of the second-order Voronoi cells listed in \code{vor2_active}
#' @field vis_par List containing the \code{par()} parameters set during a call to the the \code{vis_vor_*} functions. 
#' 
#' @section Methods:
#' Each class method has its own documentation, accessible via \code{?VorVQ::<method_name>}. 
#' For completeness, the list is repeated here in entirety.  Additional functionality 
#' for visualizing a two-dimensional VOR object is available through the \code{vis_*} functions. 
#' See their documentation for more information. 
#' 
#' \describe{
#' \item{\code{VOR$new}}{Instantiate a VOR object}
#' \item{\code{initialize_VOR}}{Initialize the VOR object with the prototypes and CADJ matrix from ANN learning}
#' \item{\code{set_bounds}}{Set the upper and lower bounds of the tessellated region}
#' \item{\code{set_vor1_active}}{Set the indices identifying the active first-order Voronoi cells}
#' \item{\code{set_vor2_active}}{Set the indices identifying the active second-order Voronoi cells}
#' \item{\code{get_params}}{Get the parameters used for various methods of the VOR object}
#' \item{\code{set_params}}{Set the parameters used for various methods of the VOR object}
#' \item{\code{calc_GADJ}}{Calculate the Gabriel adjacency matrix of the tessellation}
#' \item{\code{calc_DADJ}}{Calculate the Delaunay adjacency matrix of the tessellation}
#' \item{\code{set_DADJ}}{Set the Delaunay adjacency matrix of the tessellation, if obtained elsewhere}
#' \item{\code{get_params}}{Get the parameters used for various methods of the VOR object}
#' \item{\code{get_vor1_polytope}}{Get the polytope definition of a first-order Voronoi cell}
#' \item{\code{get_vor2_polytope}}{Get the polytope definition of a second-order Voronoi cell}
#' \item{\code{calc_vor1_centers}}{Compute the Chebhyshev centers of the active first-order Voronoi cells}
#' \item{\code{calc_vor2_centers}}{Compute the Chebhyshev centers of the active second-order Voronoi cells}
#' \item{\code{set_vor1_centers}}{Set the centers of the active first-order Voronoi cells}
#' \item{\code{set_vor2_centers}}{Set the centers of the active second-order Voronoi cells}
#' \item{\code{get_vor1_centers}}{Get the centers of the active first-order Voronoi cells}
#' \item{\code{get_vor2_centers}}{Get the centers of the active second-order Voronoi cells}
#' \item{\code{clear_vor1_centers}}{Clear the centers of the active first-order Voronoi cells}
#' \item{\code{clear_vor2_centers}}{Clear the centers of the active second-order Voronoi cells}
#' \item{\code{calc_vor1_MVIE}}{Compute the Maximum Volume Inscribed Ellipsoids (MVIE) of the active first-order Voronoi cells}
#' \item{\code{calc_vor2_MVIE}}{Compute the Maximum Volume Inscribed Ellipsoids (MVIE) of the active second-order Voronoi cells}
#' \item{\code{calc_vor1_Dikin}}{Compute the Dikin Ellipsoids of the active first-order Voronoi cells}
#' \item{\code{calc_vor2_Dikin}}{Compute the Dikin Ellipsoids of the active second-order Voronoi cells}
#' \item{\code{calc_all}}{Compute all Voronoi quantities}
#' \item{\code{save}}{Save a VOR object to disk}
#' \item{\code{load}}{Load a previously saved VOR object from disk}
#' }
#' 
NULL


# ***** Constructor *****

#' @name new
#' @title Create an empty VOR object
#' @return an empty VOR object templated for initialization (via \link{initialize_VOR}).
#' @usage VOR$new()
NULL 


# ***** Initialization *****

#' @name initialize_VOR
#' @title Initialize a VOR object
#' @description Sets up a VOR object to compute various quantities related to the Voronoi tessellation generated by a vector quantizer. 
#' @param W matrix of prototypes (codebook vectors) in rows 
#' @param CADJ the CADJ adjacency matrix produced during recall of a vector quantizer. 
#' @details This is a wrapper function to perform most steps necessary to setup a VOR object for further calculation, 
#' setting internal variables: 
#' \itemize{
#' \item \code{W}, \code{nW}
#' \item \code{d}
#' \item \code{lb}, \code{ub}
#' \item \code{CADJ}
#' \item \code{vor1_active}, \code{vor2_active}
#' \item \code{GADJ}
#' }
#' @return None
#' @usage VORobj$initialize_VOR(W, CADJ)
NULL

#' @name set_bounds
#' @title Set the bounds of the tessellated region
#' @description To ensure that all Voronoi cells are closed in \eqn{R^d} dimension-wise lower and upper 
#' bounds of the entire tessellated region must be specified. By default, the upper and lower bounds 
#' for each dimension are set = +/- 5% of \code{range(W)} in \code{initialize_VOR}.  This method is used 
#' to change these defaults to a user specified range, which must span the entire range of W 
#' in each dimension. This will be checked internally. Note that changing the global bounds can 
#' alter the Delaunay edges between prototypes which are on/near the perimiter of the point cloud defined by \code{W}. 
#' @param lb a vector of lower bounds
#' @param ub a vector of upper bounds 
#' @details Both \code{lb} and \code{ub} can have \code{length = 1}, in which case the single 
#' value will be recycled across dimension. Otherwise, they must have \code{length = ncol(W)}. 
#' Both must be supplied, even if only 
#' Internally, fields \code{lb} and \code{ub} will be overwritten. 
#' @return None
#' @usage VORobj$set_bounds(lb, ub)
NULL 

#' @name calc_GADJ
#' @title Calculate the Gabriel graph adjacency 
#' @description The Gabriel graph is a proximity graph that is a known sub-graph of the Delaunay triangulation 
#' induced by the Voronoi tessellation of the prototypes in \code{W}.  The Gabriel ADJacency matrix helps 
#' speed up the computation of the full Delaunay graph and is calculated automatically during \code{initialize_VOR} 
#' via an internal call to \code{calc_GADJ}. This method should not need to be called directly by a user but it 
#' exposed for completeness. 
#' @details Field \code{GADJ} will be overwritten internally. 
#' @return None
#' 
#' @usage VORobj$calc_GADJ()
#' 
#' @references
#' \insertRef{Gabriel1969}{VorVQ}
NULL 

#' @name get_params
#' @title Get parameters for Voronoi calculations
#' @description Several methods require parameters whose values can be accessed via this method. 
#' See \link[VorVQ]{set_params} for a description of each. 
#' @return A list with named components equal to the parameters.
#' 
#' @usage VORobj$get_params()
NULL 

#' @name set_params
#' @title Set parameters for Voronoi calculations
#' @description Several methods require parameters which can be accessed via this method. 
#' The parameters and their defaults (set during a call to \code{initialize_VOR} are: 
#' \describe{
#' \item{\code{parallel}}{boolean, whether to perform all calculations in parallel (via RcppParallel). 
#' Default = T.}
#' \item{\code{MVIE_maxiter}}{The maximum number of iterations allowed when calculating the MVIEs for Voronoi cells. Default = 100.}
#' \item{\code{MVIE_tol}}{The relative convergence tolerance for the MVIE routine, which exits when the volume of the computed MVIE changes < this tol. Default = 1e-4.}
#' \item{\code{MVIE_fix_c0}}{boolean, whether to allow the MVIE routine to update the MVIE center along with the ellipsoidal rotation. Default = T.}
#' }
#' @param param_list a list with component names in the parameter set defined above, and corresponding desired values. 
#' @return None
#' 
#' @usage VORobj$set_params(list(param_name = param_value))
NULL 


#' @name set_vor1_active
#' @title Set the active first-order Voronoi cells
#' @description The computations with a VOR object (computing Voronoi centers and ellipsoids) are only performed for "active" cells, 
#' which are identified by the prototype index (row index in \code{W}) which generated the cell. 
#' By default the active first-order Voronoi cells are those which contain data mapped to them during a recall 
#' of the ANN, as indicated by CADJ (the first-order cells corresponding to the rows of CADJ with with non-zero rowsum). 
#' A different set of active cells can be set by calling this method
#' @param vor1_indices a vector containing the indices (rows of \code{W}, 1-based) of the desired active first-order cells 
#' @return None, the field \code{vor1_active} is set internally. 
#' @usage VORobj$set_vor1_active(vor1_indices)
NULL 

#' @name set_vor2_active
#' @title Set the active second-order Voronoi cells
#' @description The computations with a VOR object (computing Voronoi centers and ellipsoids) are only performed for "active" cells, 
#' which are identified by the prototype indices (row index in \code{W}) which generated the cell. 
#' By default the active second-order Voronoi cells are those which contain data mapped to them during a recall 
#' of the ANN, as indicated by CADJ (those second-order cells corresponding to \code{CADJ(i,j) > 0}). 
#' A different set of active cells can be set by calling this method
#' @param vor2_indices a matrix whose rows contain the \code{(i,j)} indices (rows of \code{W}, 1-based) of the desired active second-order cells 
#' @return None, the field \code{vor2_active} is set internally. 
#' @usage VORobj$set_vor2_active(vor2_indices)
NULL 

#' @name calc_DADJ
#' @title Calculate the Delaunay graph adjacency 
#' @description The Delaunay graph dual of a first-order Voronoi tessellation has edges between 
#' prototype vertices whose corresponding Voronoi cells intersect (share a face). This method 
#' computes the (binary) adjacency matrix of this graph using the method Agrell (referenced below),  
#' which involves solving a linear program for each unique pair of prototypes in the tessellation. 
#' To speed up the computation, the CADJ and Gabriel adjacency matrices (set during \code{initialize_VOR}) are 
#' used to "seed" the Delaunay adjacency (as proper sub-graphs of the Delaunay graph, the existence of an edge between 
#' prototypes \code{i} and \code{j} in either necessarily indicates a corresponding edge in DADJ).  
#' Note: a Delaunay adjacency is required for all ellipsoidal calculations in VorVQ. 
#' As the many thousand LPs required for full DADJ calculation can be quite large it is recommended 
#' to perform this calculation in parallel, which can be set in \link[VorVQ]{set_params} (default is \code{parallel = TRUE}).
#' @return None, the field \code{DADJ} is set internally. 
#' @usage VORobj$calc_DADJ
#' 
#' @references
#' \insertRef{Agrell1993}{VorVQ}
NULL 

#' @name set_DADJ
#' @title Set the Delaunay graph adjacency 
#' @description 
#' A Delaunay graph adjacency is required for all ellipsoidal calculations in VorVQ. 
#' If the Delaunay adjacency is available from an external source (or calculation), it can 
#' be set with this method. Otherwise, one can be computed (and set) via \link[VorVQ]{calc_DADJ}. 
#' @return None, the field \code{DADJ} is set internally. 
#' @usage VORobj$set_DADJ
NULL 

#' @name get_vor1_polytope
#' @title Get the definition of a first-order Voronoi cell 
#' @description 
#' The definition of a first-order Voronoi cell naturally induces a half-plane representation 
#' of each cell, of the form \eqn{Ax <= b}. This method will return this half-plane definition 
#' for a particular cell. 
#' By convention, the dimension-wise lower and upper bounds stored in \code{lb} and \code{ub} 
#' are appended to the system to ensure the resulting polytope is closed. 
#' Since constraints arising from any non-Delaunay adjacent prototypes (as defined in \code{DADJ}) are redundant in the above 
#' system, these are not included in the returned polytope definition (making it non-redundant, as small as possible). 
#' @param iidx an integer identifying which first-order Voronoi cell's definition is returned. 
#' It should correspond to the row index of the generator in \code{W}.  
#' @return A list with components: 
#' \describe{
#' \item{A}{The LHS coefficient matrix of the linear inequality system}
#' \item{b}{The RHS upper bounds of the linear system, length = \code{nrow(A)}}
#' \item{cid}{A vector of indices (length = \code{nrow(A)}) indicating which prototype (generator) 
#' induced the half-plane constraint stored in the corresponding row of \code{A} and \code{b}. 
#' Constraints involving the global bounds have \code{cid = nrow(W)+1} by convention.}  
#' }
#' 
#' @usage VORobj$get_vor1_polytope(iidx)
#' 
#' @references
#' \insertRef{Agrell1993}{VorVQ}
NULL 


#' @name get_vor2_polytope
#' @title Get the definition of a second-order Voronoi cell 
#' @description 
#' The definition of a second-order Voronoi cell naturally induces a half-plane representation 
#' of each cell, of the form \eqn{Ax <= b}. This method will return this half-plane definition 
#' for a particular cell. 
#' By convention, the dimension-wise lower and upper bounds stored in \code{lb} and \code{ub} 
#' are appended to the system to ensure the resulting polytope is closed. 
#' Since constraints arising from any non-Delaunay adjacent prototypes (as defined in \code{DADJ}) are redundant in the above 
#' system, these are not included in the returned polytope definition (making it non-redundant, as small as possible). 
#' @param iidx 
#' @param jidx 
#' @details \code{iidx} and \code{jidx} are integer indices (to the generators specified in the rows of \code{W}) 
#' which, together, define the second-order Voronoi cell \code{iidx-jidx}.  
#' @return A list with components: 
#' \describe{
#' \item{A}{The LHS coefficient matrix of the linear inequality system} 
#' \item{b}{The RHS upper bounds of the linear system, length = \code{nrow(A)}} 
#' \item{cid}{A vector of indices (length = \code{nrow(A)}) indicating which prototype (generator) 
#' induced the half-plane constraint stored in the corresponding row of \code{A} and \code{b}. 
#' Constraints involving the global bounds have \code{cid = nrow(W)+1} by convention.} 
#' }
#' 
#' @usage VORobj$get_vor2_polytope(iidx, jidx)
#' 
#' @references
#' \insertRef{Agrell1993}{VorVQ}
NULL 


#' @name calc_vor1_centers
#' @title Calculate the Chebyshev centers of first-order Voronoi cells 
#' @description 
#' The \href{https://en.wikipedia.org/wiki/Chebyshev_center}{Chebyshev center} of the polytope definition 
#' (from \code{get_vor1_polytope}) for each **active** first-order Voronoi cell (as specified in \code{vor1_active}) is computed. 
#' If set, these centers are used as starting interior points for MVIE and Dikin ellipsoid calculations. 
#' Clear any previously computed centers via \code{clear_vor1_centers}. 
#' @return None, the field \code{vor1_centers} is set internally.  
#' @usage VORobj$calc_vor1_centers
NULL 

#' @name calc_vor2_centers
#' @title Calculate the Chebyshev centers of second-order Voronoi cells 
#' @description 
#' The \href{https://en.wikipedia.org/wiki/Chebyshev_center}{Chebyshev center} of the polytope definition 
#' (from \code{get_vor2_polytope}) for each **active** second-order Voronoi cell (as specified in \code{vor2_active}) is computed. 
#' These centers are used as starting interior points for MVIE and Dikin ellipsoid calculations. 
#' Clear any previously computed centers via \code{clear_vor2_centers}. 
#' @return None, the field \code{vor2_centers} is set internally.  
#' @usage VORobj$calc_vor2_centers
NULL 

#' @name set_vor1_centers
#' @title Set the centers of first-order Voronoi cells 
#' @description 
#' The MVIE and Dikin methods require a starting interior point for each active first-order Voronoi polytope. 
#' Be default, the prototypes (Voronoi cell generators, as stored in \code{W}) of each active first-order cell is used. 
#' Different interior points can be set with this method. 
#' Note: the points passed to this method **must be interior** to each first-order cell; a check will 
#' be performed and an error returned if this is not the case. 
#' Clear any previously set centers via \code{clear_vor1_centers}. 
#' @param C a matrix whose rows specify the active first-order Voronoi cell centers. Must have \code{nrow(C) = length(vor1_active)}. 
#' 
#' @return None, the field \code{vor1_centers} is set internally.  
#' @usage VORobj$set_vor1_centers
NULL 

#' @name set_vor2_centers
#' @title Set the centers of second-order Voronoi cells 
#' @description 
#' The MVIE and Dikin methods require a starting interior point for each active second-order Voronoi polytope, which 
#' can be set with this method. 
#' Note: the points passed to this method **must be interior** to each second-order cell; a check will 
#' be performed and an error returned if this is not the case. 
#' Clear any previously set centers via \code{clear_vor2_centers}. 
#' @param C a matrix whose rows specify the active second-order Voronoi cell centers. Must have \code{nrow(C) = nrow(vor2_active)}. 
#' 
#' @return None, the field \code{vor2_centers} is set internally.  
#' @usage VORobj$set_vor2_centers
NULL 

#' @name get_vor1_centers
#' @title Get the centers of first-order Voronoi cells 
#' @description 
#' If \code{calc_vor1_centers} or \code{set_vor1_centers} has been set, this method returns the centers 
#' calculated or set by these methods. Otherwise, the \code{vor1_active} prototypes in \code{W} are returned. 
#' @usage VORobj$get_vor1_centers
NULL 

#' @name get_vor2_centers
#' @title Get the centers of second-order Voronoi cells 
#' @description 
#' Returns the second-order centers that were set when calling either \code{calc_vor2_centers} or \code{set_vor2_centers}. 
#' @usage VORobj$get_vor2_centers
NULL 


#' @name clear_vor1_centers
#' @title Clear the centers of first-order Voronoi cells 
#' @description 
#' This method clears any values stored in the field \code{vor1_centers} that were previously populated by 
#' either \code{calc_vor1_centers} or \code{set_vor1_centers}. 
#' @return None, the field \code{vor1_centers} is reset to a matrix with \code{nrows = ncols = 0}. 
#' @usage VORobj$clear_vor1_centers
NULL 

#' @name clear_vor2_centers
#' @title Clear the centers of second-order Voronoi cells 
#' @description 
#' This method clears any values stored in the field \code{vor2_centers} that were previously populated by 
#' either \code{calc_vor2_centers} or \code{set_vor2_centers}. 
#' @return None, the field \code{vor2_centers} is reset to a matrix with \code{nrows = ncols = 0}. 
#' @usage VORobj$clear_vor2_centers
NULL 


#' @name calc_vor1_MVIE
#' @title Compute the MVIE for first-order Voronoi cells 
#' @description 
#' The **M**aximum **V**olume **I**nscribed **E**llipsoid inside the polytope defined by each 
#' active first-order Voronoi cell is computed and stored, according to the method of Zhang & Gao. 
#' The starting points for the routine are taken from \code{vor1_centers}, if it has been set 
#' (via either \code{calc_vor1_centers} or \code{set_vor1_centers}); otherwise the prototypes stored in \code{W} 
#' are used.  The parameters controlling the routine can be viewed or changed via \code{get/set_params}. 
#' It is recommended to allow parallel computation for MVIE calculation (i.e., \code{set_params(parallel = TRUE)}.   
#' @return None, the following fields store components of the MVIE: 
#' \describe{
#' \item{\code{vor1_MVIE_c}}{matrix of ellipsoid centers, rows correspond to entries in \code{vor1_active}}
#' \item{\code{vor1_MVIE_E}}{cube of ellipsoid rotation matrices, slices correspond to entries in \code{vor1_active}}
#' \item{\code{vor1_MVIE_logdet}}{vector of the log-determinant of the ellipsoid rotations stored in the slices of \code{vor1_MVIE_E}}
#' \item{\code{vor1_MVIE_status}}{vector of the return flags from calling the MVIE routinefor each \code{vor1_active} cell}
#' \item{\code{vor1_MVIE_volratio}}{vector of the proportional volume of each \code{vor1_active} MVIE (should sum to unity)}
#' }
#' @details 
#' Details of the MVIE routine can be found in the help of \link[VorVQ]{max_vol_inscr_ell}. 
#' 
#' @usage VORobj$calc_vor1_MVIE()
#' 
#' @references 
#' \insertRef{ZhangGao2003}{VorVQ}
NULL 


#' @name calc_vor2_MVIE
#' @title Compute the MVIE for second-order Voronoi cells 
#' @description 
#' The **M**aximum **V**olume **I**nscribed **E**llipsoid inside the polytope defined by each 
#' active second-order Voronoi cell is computed and stored, according to the method of Zhang & Gao. 
#' The starting points for the routine are taken from \code{vor2_centers}, which must be populated prior 
#' to calling this method (via either \code{calc_vor2_centers} or \code{set_vor2_centers}). 
#' The parameters controlling the routine can be viewed or changed via \code{get/set_params}. 
#' It is recommended to allow parallel computation for MVIE calculation (i.e., \code{set_params(parallel = TRUE)}.   
#' @return None, the following fields store components of the MVIE: 
#' \describe{
#' \item{\code{vor2_MVIE_c}}{matrix of ellipsoid centers, rows correspond to rows of \code{vor2_active}}
#' \item{\code{vor2_MVIE_E}}{cube of ellipsoid rotation matrices, slices correspond to rows of \code{vor2_active}}
#' \item{\code{vor2_MVIE_logdet}}{vector of the log-determinant of the ellipsoid rotations stored in the slices of \code{vor2_MVIE_E}}
#' \item{\code{vor2_MVIE_status}}{vector of the return flags from calling the MVIE routinefor each \code{vor2_active} cell}
#' \item{\code{vor2_MVIE_volratio}}{vector of the proportional volume of each \code{vor2_active} MVIE (should sum to unity)}
#' }
#' @details 
#' Details of the MVIE routine can be found in the help of \link[VorVQ]{max_vol_inscr_ell}. 
#' 
#' @usage VORobj$calc_vor2_MVIE()
#' 
#' @references 
#' \insertRef{ZhangGao2003}{VorVQ}
NULL 


#' @name calc_vor1_Dikin 
#' @title Compute the Dikin ellipsoid for first-order Voronoi cells 
#' @description 
#' The Dikin ellipsoid is a locally inscribed geometric approximator of a polytope.  
#' These ellipsoidal approximators are computed for each active first-order Voronoi cell 
#' at the points specified in \code{vor1_centers} (if these have been set); 
#' otherwise the active prototypes in \code{W} are used. 
#' @return None, the field \code{vor1_MVIE_E}, which is a cube whose slices contain the \code{vor1_active} Dikin ellipsoid 
#' rotations, is stored internally. 
#' @details 
#' Details of the Dikin ellipsoid can be found in the help of \link[VorVQ]{Dikin_Ellipsoid}. 
#' 
#' @usage VORobj$calc_vor1_Dikin()
#' 
#' @references 
#' \insertRef{Dikin1967}{VorVQ}
#' \insertRef{Boyd2004}{VorVQ}
NULL 

#' @name calc_vor2_Dikin 
#' @title Compute the Dikin ellipsoid for second-order Voronoi cells 
#' @description 
#' The Dikin ellipsoid is a locally inscribed geometric approximator of a polytope.  
#' These ellipsoidal approximators are computed for each active second-order Voronoi cell 
#' at the points specified in \code{vor2_centers} 
#' @return None, the field \code{vor2_MVIE_E}, which is a cube whose slices contain the \code{vor2_active} Dikin ellipsoid 
#' rotations, is stored internally. 
#' @details 
#' Details of the Dikin ellipsoid can be found in the help of \link[VorVQ]{Dikin_Ellipsoid}. 
#' 
#' @usage VORobj$calc_vor2_Dikin()
#' 
#' @references 
#' \insertRef{Dikin1967}{VorVQ}
#' \insertRef{Boyd2004}{VorVQ}
NULL 

#' @name calc_all
#' @title Calc all Voronoi quantities 
#' @description This wrapper method computes all Voronoi quantities exposed by the VOR class, in the following order: 
#' \itemize{
#' \item{\code{calc_DADJ}}
#' \item{\code{calc_vor2_centers}}
#' \item{\code{calc_vor1_Dikin}}
#' \item{\code{calc_vor2_Dikin}}
#' \item{\code{calc_vor1_MVIE}}
#' \item{\code{calc_vor2_MVIE}}
#' }
#' @usage VORobj$calc_all()
NULL


# ***** Saving / Loading ***** 

#' @name save 
#' @title Save a VOR object 
#' @description All fields in a VOR object can be saved to disk with this function, which allows them to be re-loaded 
#' into a new R environment at a later time (for analysis, or possibly extended training).  
#' @param vorfile a string indicating the file path and name in which to save the VOR object.  
#' This must end in extension ".vor", otherwise an error is returned. 
#' @details The VOR object is saved to disk as an R list, with each field occupying a corresponding field of the list.  
#' The file is saved in .rds format (can check its details with \link[base]{infoRDS}).  
#' 
#' Saved VORs can be re-loaded with \link{load}
#' 
#' @return None, the VOR object is saved to disk 
#' @usage VORobj$save(vorfile)
NULL 

#' @name load
#' @title Load an existing VOR object 
#' @description VOR objects previously written to disk via the \link{save} method can be re-loaded into 
#' a new R environment with this function. All fields of the internal C++ class will be populated, and 
#' all methods can be called on the loaded VOR object.  
#' @param vorfile a string indicating the file path and name of the saved VOR object.  
#' @details Because the .vor file is in .rds format it can, technically, be loaded directly into an R environment 
#' as a list via \link[base]{readRDS}. This can be useful for spot checking the contents of a saved VOR object, but 
#' does not allow use of any of its methods (or visualizations).  The \code{load} methods allows for 
#' proper restoration of a previously saved VOR 
#' 
#' @return None, the VOR object is loaded 
#' @usage VORobj$load(vorfile)
NULL 












