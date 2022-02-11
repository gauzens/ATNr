
test_that("The two versions of Scaled converge", {
  set.seed(123)
  n_species <- 50
  n_basal <- 20

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
  model <- create_model_Scaled(n_species, n_basal, masses, fw)
  # model2 uses same parameters as model
  model2 <- new(ATNr:::Scaled_loops, n_species, n_basal)
  model2[["BM"]] <- masses
  # model2[["log_BM"]] <- log10(masses)
  model2[["fw"]] <- fw
  
  biomasses <- runif(n_species, 2, 3)
  
  # 5) define the desired integration time.
  model$ext = 1e-6
  model <- initialise_default_Scaled(model)
  model$q <- 1.2

  
  # initialise properly model2
  model2$B0 = model$B0
  model2$K = model$K
  model2$X = model$X
  model2$c = model$c
  model2$e = model$e
  model2$ext = model$ext
  model2$q = model$q
  model2$r = model$r
  model2$w = model$w
  model2$max_feed = model$max_feed
  model2$alpha = model$alpha
  
  model$initialisations()
  times <- seq(0, 150000, 1500)

  sol <- lsoda_wrapper(times, biomasses, model)
  sol2 <-   deSolve::lsoda( biomasses, times,
                            func = function(t, y, pars) list(pars$ODE(y, t)),
                            model2,
  )

  extinct = tail(sol,1) < model$ext

  expect_equal(sol[nrow(sol), !extinct], sol2[nrow(sol), !extinct], tolerance = 0.0001)
})

# comparing the execution time of the two versions: 
# note: the comparison is not entirely fair. I made some optimisations in the 
# armadillo approach which is not in the loop approach

# library(rbenchmark)
# benchmark("arma: " = {sol <- lsoda_wrapper(times, biomasses, model)},
#           "loops:  " = {sol2 <- lsoda_wrapper(times, biomasses, model2)},
#           replications = 100)



