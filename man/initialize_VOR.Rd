% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vorobj_doc.R
\name{initialize_VOR}
\alias{initialize_VOR}
\title{Initialize a VOR object}
\usage{
VORobj$initialize_VOR(W, CADJ)
}
\arguments{
\item{W}{matrix of prototypes (codebook vectors) in rows}

\item{CADJ}{the CADJ adjacency matrix produced during recall of a vector quantizer.}
}
\value{
None
}
\description{
Sets up a VOR object to compute various quantities related to the Voronoi tessellation generated by a vector quantizer.
}
\details{
This is a wrapper function to perform most steps necessary to setup a VOR object for further calculation, 
setting internal variables: 
\itemize{
\item \code{W}, \code{nW}
\item \code{d}
\item \code{lb}, \code{ub}
\item \code{CADJ}
\item \code{vor1_active}, \code{vor2_active}
\item \code{GADJ}
}
}
