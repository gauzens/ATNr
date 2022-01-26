#' Default parameters as in Schneider et al. (2016)
#'
#' A dataset containing the default parameters as in the Schneider et al. (2016)
#' and used to parametrize the default models. See also
#' \code{create_model_Unscaled_nuts}, \code{create_Lmatrix},
#' \code{initialise_default_Unscaled_nuts}.
#'
#' @format A list with the default parameters:
#' \describe{ 
#' \item{Temperature}{ambient temperature in Celsius}
#' \item{T.K}{default temperature, 20 degree Celsius in Kelvin}
#' \item{k}{Boltzmann's constant}
#' \item{T0}{20 degree Celsius in Kelvin, used to estimate scaling law of metabolic rates}
#' \item{q}{Hill's exponent of the functional response}
#' \item{Ropt}{consumer/resource optimal body mass ratio}
#' \item{gamma}{shape of the Ricker function}
#' \item{mu_c}{average predator interference}
#' \item{sd_c}{standard deviation of predator interference}
#' \item{E.c}{Activation energy for interference}
#' \item{h0}{scaling constant of the power-law of handling time with consumer and resource body mass}
#' \item{hpred}{exponent associated to predator body mass for the allometric scaling of handling time}
#' \item{hprey}{exponent associated to prey body mass for the allometric scaling of handling time}
#' \item{E.h}{Activation energy for handling time}
#' \item{b0}{normalisation constant for capture coefficient}
#' \item{bprey}{exponent associated to prey body mass for the allometric scaling of capture coefficient}
#' \item{bpred}{exponent associated to predator body mass for the allometric scaling of capture coefficient}
#' \item{E.b}{Activation energy for capture coefficient}
#' \item{e_P}{Assimilation efficiency associated to the consumption of a plant species}
#' \item{e_A}{Assimilation efficiency associated to the consumption of an animal species}
#' \item{x_P}{scaling constant of the power-law of metabolic demand per unit of plant biomass}
#' \item{x_A}{scaling constant of the power-law of metabolic demand per unit of animal biomass}
#' \item{E.x}{Activation energy for metabolic rates}
#' \item{expX}{TBD}
#' \item{D}{turnover rate of nutrients}
#' \item{nut_up_min}{Minimum uptake efficiency of plants}
#' \item{nut_up_max}{Maximum uptake efficiency of plants}
#' \item{mu_nut}{Average maximum nutrient concentration}
#' \item{sd_nut}{standard deviation of maximum nutrient concentration}
#' \item{v}{relative content of nutrient 1 in plant biomass} 
#' }
#' 
#' @references Schneider, F. D., Brose, U., Rall, B. C., & Guill, C. (2016).
#'   Animal diversity and ecosystem functioning in dynamic food webs. Nature
#'   Communications, 7(1), 1-8.
"schneider"