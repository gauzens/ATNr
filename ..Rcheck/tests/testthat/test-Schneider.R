test_that("The two versions of Unscaled_nuts converge", {
  set.seed(123)
  temperature <- 20
  n_species <- 50
  n_basal <- 20
  n_nut <- 2
  #body mass of species
  masses <- 10 ^ c(sort(runif(n_basal, 1, 6)),
                   sort(runif(n_species - n_basal, 2, 9)))
  # 2) create the food web
  # create the L matrix
  
  L <- create_Lmatrix(masses,n_basal, Ropt = 100, gamma = 2, th = 0.01)
  # create the 0/1 version of the food web
  fw <- L
  fw[fw > 0] <- 1
  
  # 3) create the models:
  model <- create_model_Unscaled_nuts(n_species, n_basal, n_nut, masses, fw)
  # model2 uses same parameters as model
  model2 <- new(Unscaled_nuts_loops, n_species, n_basal, n_nut)
  model2[["BM"]] <- masses
  model2[["log_BM"]] <- log10(masses)
  model2[["fw"]] <- fw
  
  biomasses <- runif(n_species + n_nut, 2, 3)
  
  # 5) define the desired integration time.
  model$ext = 1e-6
  model <- initialise_default_Unscaled_nuts(model, L, temperature = temperature)
  model$q <- 0.2
  model$S <- rep(60, n_nut)
  
  # initialise properly model2
  model2$D = model$D
  model2$K = model$K
  model2$S = model$S
  model2$X = model$X
  model2$V = model$V
  model2$b = model$b
  model2$c = model$c
  model2$e = model$e
  model2$ext = model$ext
  model2$h = model$h
  model2$q = model$q
  model2$r = model$r
  model2$w = model$w
  
  model$initialisations()
  times <- seq(0, 150000, 1500)
  sol <- lsoda_wrapper(times, biomasses, model)
  sol2 <- lsoda_wrapper(times, biomasses, model2)
  
  extinct <- tail(sol, 1) < model$ext
  
  expect_equal(sol[nrow(sol), !extinct], sol2[nrow(sol), !extinct], tolerance = 0.0001)
})

# comparing the execution time of the two versions: 
# note: the comparison is not entirely fair. I made some optimisation in the 
# armadillo approach which is not in the loop approach
# library(rbenchmark)

# benchmark("arma: " = {sol <- lsoda_wrapper(times, biomasses, model)},
#           "loops:  " = {sol2 <- lsoda_wrapper(times, biomasses, model2)},
#           replications = 100)
