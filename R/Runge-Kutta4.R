#' @keywords internal
.runge_kutta4<- function(t, biomasses, model){
  bioms <- matrix(NA, ncol = length(biomasses), nrow = length(t))
  biom.step<- biomasses
  delta.t <- t[2] - t[1]
  for (i in 1:length(t)){
    bioms[i, ] <- biom.step
    k1 <- model$ODE(biom.step, i*delta.t)
    k2 <- model$ODE(biom.step + 0.5*delta.t * k1, (i+0.5)*delta.t)
    k3 <- model$ODE(biom.step + 0.5*delta.t * k2, (i+0.5)*delta.t)
    k4 <- model$ODE(biom.step + delta.t * k3, i*delta.t)
    biom.step <- biom.step + (delta.t/6) * (k1 + 2*k2 + 2*k3 + k4)
  }
  return(cbind(t, bioms))
}
