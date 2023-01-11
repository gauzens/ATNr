#' @title Plot food web dynamics
#'
#' @description Plot solution of the ODE for the food web. Currently only
#'   species and not nutrients are plotted.
#'
#' @param x matrix with solutions. First row should be the time vector.
#' @param nb_s numeric, number of species as in the model (e.g.,
#'   \code{create_model_Unscaled_nuts}).
#'   
#' @return No return value, called for side effects.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ATNr)
#' library(deSolve)
#' set.seed(123)
#' # number of species, nutrients, and body masses
#' n_species <- 20
#' n_basal <- 5
#' n_nutrients <- 3
#' masses <- sort(10^runif(n_species, 2, 6)) #body mass of species
#' # create food web matrix
#' L <- create_Lmatrix(masses, n_basal)
#' L[, 1:n_basal] <- 0
#' fw <- L
#' fw[fw > 0] <- 1
#' model <- create_model_Unscaled_nuts(
#'   n_species,
#'   n_basal,
#'   n_nutrients,
#'   masses,
#'   fw
#' )
#' # initialize model as default in Schneider et al. (2016)
#' model <- initialise_default_Unscaled_nuts(model, L)
#' # defining integration time
#' times <- seq(0, 500, 5)
#' biomasses <- runif(n_species + n_nutrients, 2, 3)
#' sol <- lsoda_wrapper(times, biomasses, model, verbose = FALSE)
#' plot_odeweb(sol, model$nb_s)
#' }
plot_odeweb <- function(x, nb_s) {
  stopifnot((ncol(x) - 1) >= nb_s)
  pal <- grDevices::colorRampPalette(c("blue", "red"))(nb_s)
  pal <- grDevices::adjustcolor(pal, alpha.f = .5)
  plot(c(0, max(x[, 1])), #xlim
       c(0, max(x[, c((ncol(x) - nb_s + 1) : ncol(x))])), #ylim
       frame = FALSE,
       xlab = "Time",
       ylab = "Biomass",
       col = NA)
  for (i in seq(ncol(x) - nb_s + 1, ncol(x))) {
    graphics::points(x[, 1], x[, i], col = pal[i - ncol(x) + nb_s], pch = 20, cex = .5)
    graphics::lines(x[, 1], x[, i], col = pal[i - ncol(x) + nb_s], lw = 1)
  }
}
