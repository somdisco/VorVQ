# VorVQ

`VorVQ` exposes a collection of tools for computing quantities which describe the geometry of the cells of the first and second-order Voronoi tessellations induced by the learned prototypes (i.e., codebook vectors) of a vector quantizer.  These geometric descriptions are formed via computation of:

+ Half-plane representations of the closed, convex polytopes comprising the Voronoi tessellation (i.e., the Voronoi cells, which we also call Voronoi polytopes in what follows).
+ The Chebyshev centers [@Boyd2004] of each Voronoi polytope.
+ The adjacency matrix of the Delaunay triangulation [@Delaunay1934] of the prototype vectors in $\mathbb{R}^d$, for arbitrary dimension $d$.  The well-known [`Qhull`](http://www.qhull.org/) library (available for R in the [`geometry`](https://cran.r-project.org/web/packages/geometry/) package) computed Delaunay triangulations for low ($\lessapprox 5$) dimension; to our knowledge `VorVQ` is the only package to allow equivalent computation in high dimension. + The adjacency matrix of the Gabriel graph ([@Gabriel1969]), which is a [known sub-graph](https://en.wikipedia.org/wiki/Gabriel_graph) of the Delaunay triangulation. The R package [cccd](https://cran.r-project.org/web/packages/cccd/) also exposes this functionality; it is re-cast here with support for parallel computation to expedite the calculation of the Delaunay triangulation. 
+ The Maximum Volume Inscribed ([@ZhangGao2003]) and Dikin ([@Dikin1967]) ellipsoids for each Voronoi polytope, which are useful for approximating the orientation and size (volume) of each Voronoi cell.  

Calculation within `VorVQ` is implemented in C++ (via [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)) with multi-threaded support (via [RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)). 

# Installation

```r
devtools::install_github("somdisco/VorVQ")
```

`VorVQ` requires a valid installation of the [Gurobi](https://www.gurobi.com/) optimization library. Free academic licenses are available at Gurobi's website.  


# Documentation

See the [VorVQ homepage](https://somdisco.github.io/VorVQ/output/index.html) for more information.
