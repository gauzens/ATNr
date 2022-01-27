#' @title Function to remove species from a model class
#'
#' @param species integer vector, the indices of species to remove.
#' @param model model object
#' @param nuts integer vector, the indices of nutrients to remove. Parameter
#'   specific to the Unscaled_nuts model.
#' @export
#' @return A model object where the data structure has bee updated to remove the
#'   species in parameters.
remove_species = function(species, model, nuts = NULL){
  if (class(model)[1] %in% c("Rcpp_Scaled", "Rcpp_Scaled_loops")){
    model2 = .remove_Scaled(species, model)
  } else if (class(model)[1] %in% c("Rcpp_Unscaled_nuts", "Rcpp_Unscaled_nuts_loops")) {
    model2 = .remove_Unscaled_nuts(species, model, nuts)
  } else if (class(model)[1] %in% c("Rcpp_Scaled", "Rcpp_Scaled_loops")) {
    model2 = .remove_Scaled(species, model)
  } else {
    stop("model must be a model class from the ATNr package")
  }
  return(model2)
}

#' Internal function to remove species from a Scaled model
#'
#'@keywords internal
#'
.remove_Scaled <- function(species, model) {
  # consumers: indices of consumer species in data that does not have basal species
  consumers <- species[species > model$nb_b] - model$nb_b
  # index of basal species
  basals <- species[species <= model$nb_b]

  new.nb_b <- model$nb_b - sum(species < model$nb_b)
  new.nb_s <- model$nb_s - length(species)

  model2 <- methods::new(Scaled, new.nb_s, new.nb_b)

  model2$fw <- model$fw[-species, -species]
  model2$BM <- model$BM[-species]
  model2$X <- model$X[-species]
  model2$e <- model$e[-species]
  model2$log_BM <- model$log_BM[-species]
  model2$c <- model$c
  model2$q <- model2$q

  model2$dB <- model$dB[-species]

  # mat[-x,-y] or vec[-x] return void objects when x or y are numeric(0).
  # then I guess I have no other choice than checking consumers and basals?
  if (length(consumers) > 0) {
    model2$B0 <- model$B0[-consumers]
    model2$F <- model$F[-species, -consumers]
    model2$max_feed <- model$max_feed[-consumers]
    model2$w <- model$w[-species, -consumers]
  }else{
    model2$B0 <- model$B0
    model2$F <- model$F[-species, ]
    model2$max_feed <- model$max_feed
    model2$w <- model$w[-species, ]
  }
  # same checks for basals
  if (length(basals) > 0) {
    model2$r <- model$r[-basals]
  }else{
    model2$r <- model$r
  }

  return(model2)
}


#' Internal function to remove species from a Unscaled_nuts model
#'
#'@keywords internal
.remove_Unscaled_nuts <- function(species, model, nuts) {
  # consumers: indices of consumer species in data that does not have basal species
  consumers <- species[species > model$nb_b] - model$nb_b
  # index of basal species
  basals <- species[species <= model$nb_b]

  new.nb_b <- model$nb_b - sum(species < model$nb_b)
  new.nb_s <- model$nb_s - length(species)

  new.nb_n <- model$nb_n
  if (length(nuts > 0)) {
    new.nb_n <- model$nb_n - 1
  }

  model2 <- methods::new(Unscaled_nuts, new.nb_s, new.nb_b, new.nb_n)

  model2$fw <- model$fw[-species, -species]
  model2$BM <- model$BM[-species]
  model2$X <- model$X[-species]
  model2$e <- model$e[-species]
  model2$log_BM <- model$log_BM[-species]
  model2$c <- model$c
  model2$q <- model2$q

  model2$dB <- model$dB[-(species + model$nb_n)]
  model2$dB <- model$dB[-nuts]

  # mat[-x,-y] or vec[-x] return void objects when x or y are numeric(0).
  # then I guess I have no other choice than checking consumers and basals?
  if (length(consumers) > 0){
    model2$F <- model$F[-species, -consumers]
    model2$b <- model$b[-species,-consumers]
    model2$w <- model$w[-species, -consumers]
    model2$h <- model$h[-species,-consumers]
  }else{
    model2$F <- model$F[-species, ]
    model2$w <- model$w[-species,]
    model2$b <- model$b[-species,]
    model2$h <- model$h[,-consumers]
    
  }
  # same checks for basals
  if (length(basals) > 0){
    model2$r <- model$r[, -basals]
    model2$K <- model$K[, -basals]
    model2$V <- model$V[, -basals]
    model2$uptake <- model$uptake[, -basals]
  }else{
    model2$r <- model$r
    model2$K <- model$K
    model2$V <- model$V
    model2$uptake <- model$uptake
  }
  # now nutrients, no need to check for plants, as already done before

  if (!is.null(nuts)) {
    model2$K <- model$K[-nuts, ]
    model2$S <- model$S[-nuts]
    model2$V <- model$V[-nuts, ]
  }

  return(model2)
}


#' Internal function to remove species from a Unscaled model
#'
#' @keywords internal
.remove_Unscaled = function(species, model) {
  # consumers: indices of consumer species in data that does not have basal species
  consumers <- species[species > model$nb_b] - model$nb_b
  # index of basal species
  basals <- species[species <= model$nb_b]

  new.nb_b <- model$nb_b - sum(species < model$nb_b)
  new.nb_s <- model$nb_s - length(species)

  model2 <- methods::new(Scaled, new.nb_s, new.nb_b)

  model2$fw <- model$fw[-species, -species]
  model2$BM <- model$BM[-species]
  model2$X <- model$X[-species]
  model2$e <- model$e[-species]
  model2$log_BM <- model$log_BM[-species]
  model2$c <- model$c
  model2$q <- model2$q

  model2$dB <- model$dB[-species]

  # mat[-x,-y] or vec[-x] return void objects when x or y are numeric(0).
  # then I guess I have no other choice than checking consumers and basals?
  if (length(consumers) > 0) {
    model2$F <- model$F[-species, -consumers]
    model2$a <- model$a[-species, -consumers]
    model2$h <- model$h[-species, -consumers]
  }else{
    model2$F <- model$F[-species, ]
    model2$b <- model$b[-species, ]
    model2$h <- model$h[, -consumers]
  }
  # same checks for basals
  if (length(basals) > 0) {
    model2$r <- model$r[, -basals]
    model2$K <- model$K[-basals]
  }else{
    model2$r <- model$r
    model2$K <- model$K
  }
  
  return(model2)
}
