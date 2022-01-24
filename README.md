# To do

Add the binzer model in the run_checks function

There might be some oportunities to optimise execution time of the armadillo models, especially regarding bracketing of expressions based on matrix algebra. Not important to investigate, but could be nice to have a look at some point. Need a good knowledge of the approach though

Expend the different examples in the vignette: do 1 ecological example per model: Schneider: temperature, Delmas: complexity, Binzer: pred-prey body mass ratios?

As the models coded here differ a bit from the different published version (because some approximations where corrected or because we aim for something more generic), we need to make a document presenting the equations of the different models. A latex vignette that can be released with the package maybe? 

Incorporate the generalisation of carrying capacity presented in Delmas' paper (to select whether or not plant share the same nutrient pool) in the Delmas and Binzer classes 

Some parameters now use the same values for all nodes it refers to. Maybe we should change that to vectors to make them node specific and having something more general: 
* In schneider: D (turnover rate of the nutrients) is a single value common to all nutrient pools. 
* In all models, interspecific competition is a scalar. 


Schneider uses some normal distribution to generate parameter values. We should move to truncated normal distributions to ensure that all parameters are > 0

For the model based on Armadillo, I made a function "initialisation" that pre-compute variables needed in the ODE calculation but that do not change over time. FOr instance, the mulitplication of attack rate and handling time does not need to be calculated at each time step, as it is independant of species biomasses. Therefore, before launching an integration, the initialisation funciton has to be called. This is made to improve calculation time, but can be very error prone from an user point of view. We have to decide whether we want to keep this approach or to put all of these calculations in the ODE function (at the cost of calculation time). 

In the vignette describing the model, associate the mathematical variables to the names used in the code (in the tables)

Document the Binzer class

Add the references for the units in the model description vignette

rename the package to ATNr?
