#' @title Run checks on model parameters
#'
#' @description Check if the dimensions of vectors and matrices used in the model are correct.
#' If any dimension is not correct, an error message is returned.
#'
#' @param model a model object.
#' @param verbose Boolean, whether a message should be printed when all checks were successful
#'
#' @return No return value, only throw an error if parameters are inconsistent. 
#' 
#' @export
run_checks <- function(model, verbose = TRUE) {
  if (class(model)[1] == "Rcpp_Scaled") {
    with(model, {
      if (length(X) != nb_s) stop(" vector of metabolic rates ($X) mispecified")
      else if (length(max_feed) != nb_s - nb_b) stop(" vector of maximum feeding rates ($max_feed) mispecified")
      else if (length(e) != nb_s) stop(" vector of assimilation efficiencies mispecified")
      else if (length(r) != nb_b) stop(" vector of plant growth rate mispecified")
      else if (length(BM) != nb_s) stop(" vector of body masses ($BM) mispecified")
      else if (length(B0) != nb_s - nb_b) stop(" vector of half saturation densities ($B0) mispecified")
      else if (any(dim(fw) != c(nb_s, nb_s))) stop(" food web ($fw) dimension is incorrect")
      else if (any(dim(w) != c(nb_s, nb_s - nb_b))) stop(" dimensions of $w are incorrect")
      else if (any(dim(alpha) != c(nb_b, nb_b))) stop(" dimensions of plant competition matrix ($alpha) are incorrect")
    })
    if (verbose) message("All checks successful")
  }
  else if (class(model)[1] == "Rcpp_Unscaled_nuts") {
    with(model, {
      if (any(dim(K) != c(nb_n, nb_b))) stop(" vector of plant half saturation densities ($K2) mispecified")
      else if (length(e) != nb_s) stop(" vector of assimilation efficiencies ($e) mispecified")
      else if (length(r) != nb_b) stop(" vector of plant growth rate ($r) mispecified")
      else if (length(BM) != nb_s) stop(" vector of body masses ($BM) mispecified")
      else if (length(S) != nb_n) stop(" vector of maximal nutrient level ($S) mispecified")
      else if (any(dim(fw) != c(nb_s, nb_s))) stop(" food web dimension ($fw) is incorrect")
      else if (any(dim(w) != c(nb_s, nb_s - nb_b))) stop(" diemsions of w are incorrect")
      else if (any(dim(b) != c(nb_s, nb_s - nb_b))) stop(" dimensions of atack rates matrix ($b) are incorrect")
      else if (any(dim(h) != c(nb_s, nb_s - nb_b))) stop(" dimensions of handling times matrix ($h) are incorrect")
      else if (length(X) != nb_s) stop(" vector of metabolic rates ($X) mispecified")
      else if (any(dim(V) != c(nb_n, nb_b))) stop(" matrix relative content in the plant species' biomass ($V) mispecified")
    })
    if (verbose) message("All checks successfull")
  }
  else if (class(model)[1] == "Rcpp_Unscaled") {
    with(model, {
      if (length(e) != nb_s) stop(" vector of assimilation efficiencies ($e) mispecified")
      else if (length(r) != nb_b) stop(" vector of plant growth rate ($r) mispecified")
      else if (length(BM) != nb_s) stop(" vector of body masses ($BM) mispecified")
      else if (any(dim(fw) != c(nb_s, nb_s))) stop(" food web dimension ($fw) is incorrect")
      else if (any(dim(a) != c(nb_s, nb_s - nb_b))) stop(" dimensions of atack rates matrix ($a) are incorrect")
      else if (any(dim(h) != c(nb_s, nb_s - nb_b))) stop(" dimensions of handling times matrix ($h) are incorrect")
      else if (length(X) != nb_s) stop(" vector of metabolic rates ($X) mispecified")
      else if (any(dim(alpha) != c(nb_b, nb_b))) stop(" dimensions of plant competition matrix ($alpha) are incorrect")
    })
    if (verbose) message("All checks successfull")
  }
  # The following models are present for the testing unit only. NOt to be used
  else if (!class(model) %in% c("Rcpp_Scaled_loops", "Rcpp_Unscaled_nuts_loops", "Rcpp_Unscaled_loops") ) stop(class(model)[1], " is not supported.")
}
