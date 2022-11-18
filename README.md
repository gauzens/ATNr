Beta version of the ATNr package that propose solutions to estimate populations dynamics in food webs with Allometric Trophic Networks

# To do

There might be some opportunities to optimise execution time of the armadillo models, especially regarding bracketing of expressions based on matrix algebra. Not important to investigate, but could be nice to have a look at some point. Need a good knowledge of the approach though

Some parameters now use the same values for all nodes it refers to. Maybe we should change that to vectors to make them node specific and having something more general: 
* In Schneider: D (turnover rate of the nutrients) is a single value common to all nutrient pools. 
* In all models, interspecific competition is a scalar. 
 
