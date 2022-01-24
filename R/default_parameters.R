# Create default parameters data.
# This is only for development purporses: do not use directly.

###############################################################
##### definition of global parameters #########################
###############################################################

# Temperature <- 20
# schneider <- list(
#   Temperature = Temperature,
#   T.K = Temperature + 273.15,
#   k = 8.6173324e-5,
#   T0 = 293.15,
#   
#   q = 1.2,
#   
#   # used to generate L, foodWeb, bm
#   # optimal consumer-resource body-mass ratio
#   Ropt = 100,
#   # Ricker function width
#   gamma = 2,
#   
#   # used in feeding function
#   # consumer interference, proportion of time a consumer spences encountering con-specifics
#   mu_c = 0.8,
#   sd_c = 0.2,
#   E.c = -0.65,
#   # used in Handling Time (h)
#   h0 = 0.4,
#   hpred = rnorm(1,-0.48,0.03),
#   hprey = rnorm(1,-0.66,0.02),
#   E.h = 0.26,
#   
#   # feeding rates:
#   b0 = 50,
#   bprey = rnorm(1, 0.15, 0.03),
#   bpred = rnorm(1, 0.47, 0.04),
#   E.b = -0.38,
#   # used in calculating change in biomass of Animal or Plant
#   # NOTE called conversion efficiency in Schneider et al., 2016
#   # assimilation efficiency herbivore
#   e_P = 0.545, # NOTE From Lang et al., 2017 ORIGINAL value 0.45 from Schneider et al., 2016
#   # assimilation efficiency carnivore
#   e_A = 0.906, # NOTE From Lang et al., 2017 ORIGINAL value  0.85 from Schneider et al., 2016
#   # metabolic rate scaling constant of plants
#   x_P = 0.141, # ORIGINAL value 0.138 from Schneider et al., 2016
#   # metabolic rate scaling constant of animals
#   x_A = 0.314,  # NOTE ORIGINAL value 0.314 from Schneider, 0.305 From Ehmes et al., 2011
#   E.x = -0.69,
#   expX = -0.305,
#   
#   #used in calculating change in nutrient concentration
#   # global nutrient turn over rate (rate of replenishment)
#   D = 0.25,
#   # min and max nutrient uptake efficiencies
#   nut_up_min = 0.1,
#   nut_up_max = 0.2,
#   # nutrient 'densities'
#   mu_nut = 10,
#   sd_nut = 2,
#   # plant nutrient proportions: 
#   v1 = 1,
#   v2 = 0.5
# )
# 
# usethis::use_data(schneider, overwrite = TRUE)
