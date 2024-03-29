% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_model.R
\name{create_model_Scaled}
\alias{create_model_Scaled}
\title{Initialize an ATN model, following Delmas et al. 2017, Methods in Ecology and Evolution}
\usage{
create_model_Scaled(nb_s, nb_b, BM, fw)
}
\arguments{
\item{nb_s}{integer, number of total species.}

\item{nb_b}{integer, number of basal species.}

\item{BM}{float vector, body mass of species.}

\item{fw}{binary adjacency matrix of the food web.}
}
\value{
An object of class \emph{ATN (Rcpp_parameters_prefs)}.
}
\description{
Initialize an ATN model, following Delmas et al. 2017, Methods in Ecology and Evolution
}
\details{
A model is defined by the total number of species
(\emph{nb_s}), the number of basal species (\emph{nb_b}),
the number of nutrients (\emph{nb_n}), the body masses
(\emph{BM}) of species, and the adjacency matrix (\emph{fw})
representing species interactions.
}
\examples{
library(ATNr)
set.seed(123)
n_species <- 50
n_basal <- 20
masses <- sort(10^runif(n_species, 2, 6)) #body mass of species
L <- create_Lmatrix(masses, n_basal)
fw <- L
fw[fw > 0] <- 1
mod <- create_model_Scaled(n_species, n_basal, masses, fw)
}
\references{
Delmas, E., Brose, U., Gravel, D., Stouffer, D.B. and Poisot, T.
  (2017), Simulations of biomass dynamics in community food webs. Methods
  Ecol Evol, 8: 881-886. https://doi.org/10.1111/2041-210X.12713
}
