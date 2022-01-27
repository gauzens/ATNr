#' @title Estimate the Jacobian matrix of a ODE system
#'
#' @param bioms float vector, biomass of species.
#' @param ODE function that computes the ODEs from one of the model available
#' @param eps float, scale precision of the numerical approximation.
#' @export
#' @details The function provides a numerical estimation of the Jacobian matrix
#' based on the 5 points stencil method. The precision of the method is in  \deqn{O(h^5)},
#' where \deqn{h = eps*bioms}. The choice of eps should ensure that \deqn{h^5}
#' is always lower to the extinction threshold.
#'
#' The dimension of the Jacobian matrix are not always matching the number of species in the system.
#' This is because we considered that a perturbation can not correspond to the recolonisation of an extinct species.
#' Therefore, extinct species are removed from the system to calculate the Jacobian matrix.
#' @return A matrix corresponding to the Jacobian  of the system estimated at the parameter biomasses

Joacobian <- function(bioms, ODE, eps = 1e-8){
  nb_s <-  length(bioms)
  Jacob <- matrix(NA, nrow = nb_s, ncol = nb_s)
  for (cons in 1:nb_s){
    # h:  magnitude of the perturbation
    h <- eps * bioms[cons];
    # generate the vectors with a slight variation applied
    hs <- bioms
    hs2 <- bioms
    mhs <- bioms
    mhs2 <- bioms
    # here I consider that extinct species can't perturb the system (they are gone)
    # so I apply he perturbation only for non extinct species
    if (bioms[cons] > 0){
      hs[cons] <- hs[cons] + h;
      hs2[cons] <- hs2[cons] + 2 * h;
      mhs[cons] <- mhs[cons] - h;
      mhs2[cons] <- mhs2[cons] - 2 * h;
    }
    # compute the local derivatives
    res.values <- (-ODE(hs2, 0.0) + 8*ODE(hs, 0.0) - 8*ODE(mhs, 0.0) + ODE(mhs2, 0.0)) / (12*h)
    # fill the jacobian matrix
    Jacob[, cons] <- res.values

  }
  return(Jacob[bioms > 0, bioms > 0])
}
