% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connected_component.R
\name{is_connected}
\alias{is_connected}
\title{Detect whether a food web is composed of several disconnected sub-networks}
\usage{
is_connected(fw)
}
\arguments{
\item{fw}{binary adjacency matrix of the food web.}
}
\value{
Boolean: TRUE if the food web is connected, FALSE if several disconnected sub-networks are detected.
}
\description{
Run a deep search first algorithm (DFS)
}
\examples{

library(ATNr)
set.seed(123)
# number of species, nutrients, and body masses
n_species <- 20
n_basal <- 5
n_nutrients <- 3
masses <- sort(10^runif(n_species, 2, 6)) #body mass of species
# create food web matrix
L <- create_Lmatrix(masses, n_basal)
L[, 1:n_basal] <- 0
fw <- L
fw[fw > 0] <- 1
connected <- is_connected(fw)
}
