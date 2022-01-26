Beta version of the ATNr package that propose solutions to estimate populations dynamics in food webs with Allometric Trophic Networks

# To do


There might be some opportunities to optimise execution time of the armadillo models, especially regarding bracketing of expressions based on matrix algebra. Not important to investigate, but could be nice to have a look at some point. Need a good knowledge of the approach though

Some parameters now use the same values for all nodes it refers to. Maybe we should change that to vectors to make them node specific and having something more general: 
* In schneider: D (turnover rate of the nutrients) is a single value common to all nutrient pools. 
* In all models, interspecific competition is a scalar. 

Schneider uses some normal distribution to generate parameter values. We should move to truncated normal distributions to ensure that all parameters are > 0

For the model based on Armadillo, I made a function "initialisation" that pre-compute variables needed in the ODE calculation but that do not change over time. For instance, the multiplication of attack rate and handling time does not need to be calculated at each time step, as it is independant of species biomasses. Therefore, before launching an integration, the initialisation function has to be called. This is made to improve calculation time, but can be very error prone from an user point of view. We have to decide whether we want to keep this approach or to put all of these calculations in the ODE function (at the cost of calculation time). 

Add the references for the units in the model description vignette
