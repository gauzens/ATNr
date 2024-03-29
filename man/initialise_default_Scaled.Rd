% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_initialisations.R
\name{initialise_default_Scaled}
\alias{initialise_default_Scaled}
\title{Default parameters for the scaled version of ATN as in Delmas et al.
  2016}
\usage{
initialise_default_Scaled(model)
}
\arguments{
\item{model}{an object of class \emph{Rcpp_Scaled}.}
}
\value{
An object of class \emph{Rcpp_Scaled} with default
  parameters as in Delmas et al. (2017).
}
\description{
Initialise the default parametrisation for the scaled version of
  the ATN model as in Delmas et al. (2016).
}
\examples{
library(ATNr)
set.seed(123)
masses <- runif(20, 10, 100) #body mass of species
L <- create_Lmatrix(masses, 10, Ropt = 10)
L[L > 0] <- 1
mod <- create_model_Scaled(20, 10, BM = masses, fw = L)
mod <- initialise_default_Scaled(mod)
}
\references{
Delmas, E., Brose, U., Gravel, D., Stouffer, D.B. and Poisot, T.
  (2017), Simulations of biomass dynamics in community food webs. Methods
  Ecol Evol, 8: 881-886. https://doi.org/10.1111/2041-210X.12713
}
