# library(ATNr)
# devtools::load_all()
# devtools::run_examples()
# devtools::build_vignettes()
devtools::test()



# 
# set.seed(123)
# n_species <- 50
# n_basal <- 20
# 
# #body mass of species
# masses <- 10 ^ c(sort(runif(n_basal, 1, 6)),
#                  sort(runif(n_species - n_basal, 2, 9)))
# # 2) create the food web
# # create the L matrix
# L <- create_Lmatrix(masses,n_basal, Ropt = 100, gamma = 2, th = 0.01)
# # create the 0/1 version of the food web
# fw <- L
# fw[fw > 0] <- 1
# 
# # 3) create the models:
# model <- create_model_Unscaled(n_species, n_basal, masses, fw)
# 
# 
# biomasses <- runif(n_species, 30, 35)
# 
# # 5) define the desired integration time.
# model$ext = 1e-6
# model <- initialise_default_Unscaled(model)
# model$q <- 0.2
# 
# model$initialisations()
# times <- seq(0, 355000000, 1000000)
# 
# sol <- lsoda_wrapper(times, biomasses, model)
# 
# 
# model <- create_model_Scaled(n_species, n_basal, masses, fw)
# 
# model$ext = 1e-6
# model <- initialise_default_Scaled(model)
# model$q <- 1.2
# 
# 
# model$initialisations()
# times <- seq(0, 150000, 1500)
# 
# sol <- lsoda_wrapper(times, biomasses, model)

