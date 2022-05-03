#' @title Wrapper for lsoda
#'
#' @description This is a wrapper to call \code{lsoda} from
#' \emph{deSolve} and solve the ODE.
#' Package \code{deSolve} needs to be installed to run
#' this wrapper.
#'
#' @param t vector of times.
#' @param y vector of biomasses.
#' @param model object of class \emph{ATN (Rcpp_parameters_prefs)}.
#' @param verbose Boolean, whether a message should be printed when all checks were successful 
#' @param ... additional arguments to pass to `lsoda`
#' @export
#'
#' @return A matrix for the ODE solution with species as columns and
#' times as rows.
#'
#' @examples
#' library(ATNr)
#' set.seed(123)
#' masses <- runif(20, 10, 100) #body mass of species
#' L <- create_Lmatrix(masses, 10, Ropt = 10)
#' L[L > 0] <- 1
#' mod <- create_model_Unscaled_nuts(20, 10, 3, masses, L)
#' mod <- initialise_default_Unscaled_nuts(mod, L)
#' biomasses <- masses ^ -0.75 * 10 ^ 4 #biomasses of species
#' biomasses <- append(runif(3, 20, 30), biomasses)
#' times <- seq(0, 100, 1)
#' sol <- lsoda_wrapper(times, biomasses, mod)
lsoda_wrapper <- function(t, y, model, verbose = FALSE, ...) {
  model$initialisations()
  run_checks(model, verbose)
  deSolve::lsoda(
    y,
    t,
    func = function(t, y, pars) list(pars$ODE(y, t)),
    model,
    ...
  )
}

# #' @title Wrapper for sundial
# #' 
# #' @keywords internal
# #'
# #' @description This is a wrapper to call \code{cvode} from
# #' \emph{sundialr} and solve the ODE.
# #' Package \code{sundialr} needs to be installed to run
# #' this wrapper.
# #'
# #' @param t vector of times.
# #' @param y vector of biomasses.
# #' @param model object of class \emph{ATN (Rcpp_parameters_prefs)}.
# #'
# #' @return A matrix for the ODE solution with species as columns and
# #' times as rows.
# #'
# #' @examples
# #' library(ATNr)
# #' masses <- runif(50, 10, 100) #body mass of species
# #' L <- create_Lmatrix(masses, 10, Ropt = 50)
# #' L[L > 0] <- 1
# #' mod <- create_model_Unscaled_nuts(20, 10, 3, masses, L)
# #' mod <- initialise_default_Unscaled_nuts(mod, L)
# #' biomasses <- masses ^ -0.75 * 10 ^ 4 #biomasses of species
# #' biomasses <- append(runif(3, 20, 30), biomasses)
# #' times <- seq(0, 100, 1)
# #' sol <- sundial_wrapper(times, biomasses, mod)
# #' t <- times
# #' y <- biomasses
# 
# #' sundialr::cvode(
#   #' time_vector = 0.0, #time vectors
#   #' IC = y, #initial conditions
#   #' input_function = function(t, y, p) mod$ODE(y, t), #anonymous function to reorder input
#   #' Parameters = c(0, 0) #this does nothing, but is necessary for compatibility
# #' )
# 
# sundial_wrapper <- function(t, y, model) {
#   wrapper.ODE <- function(t, y, p) {
#     return(model$ODE(y, t))
#   }
#   sundialr::cvode(
#     t,
#     y,
#     input_function = wrapper.ODE,
#     c(0, 0)
#   )
# }
