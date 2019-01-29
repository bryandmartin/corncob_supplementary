args <- commandArgs(TRUE)

## Parse the arguments
if (length(args) == 0) {
  print("No arguments supplied.")
} else {
  for (i in 1:length(args)) {
    eval(parse(text = args[i]))
  }
}



devtools::install_github("bryandmartin/corncob")
require(GenSA)
require(corncob)
require(VGAM)
require(knitr)
require(ggplot2)
require(magrittr)
require(dplyr)
require(doParallel)
require(doSNOW)
require(simulator)
require(brglm2)
require(Rmpi)

# This is the main simulator file

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

# B now in method_functions to make draws more efficient
nsim <- ceiling(10000/100)

mysocks <- rep(c("b09.local", "b10.local", "b11.local", "b12.local", "b13.local",
                 "b14.local", "b15.local", "b16.local", "b17.local", "b18.local"), each = 10)

if (arg1 == 1) {
  if (!file_test("-f", "draws1.rda")) {
    draws1 <-  new_simulation(name = "setting1",
                              label = "Setting 1") %>%
      generate_model(make_betabin, seed = 1,
                     n = list(10,30,100), 
                     setting = 1,
                     vary_along = "n") %>%
      simulate_from_model(nsim = nsim, index = 1:100)
    save(draws1, file = "draws1.rda")
  }
  load("draws1.rda")
  sim1 <- draws1 %>%
    run_method(setting1, parallel = list(socket_names = mysocks,
                                         libraries = "corncob")) %>%
    evaluate(list(n, wald, pbwald, lrt, pblrt))
  save(sim1, file = paste("TypeI_setting1_", nsim*100, "_", B, ".rda",
                          sep = ""))
} else if (arg1 == 2) {
  if (!file_test("-f", "draws2.rda")) {
    draws2 <- new_simulation(name = "setting2",
                             label = "Setting 2") %>%
      generate_model(make_betabin, seed = 1,
                     n = list(10,30,100), 
                     setting = 2,
                     vary_along = "n") %>%
      simulate_from_model(nsim = nsim, index = 1:100)
    save(draws2, file = "draws2.rda")
  }
  load("draws2.rda")
  sim2 <- draws2 %>%
    run_method(setting2, parallel = list(socket_names = mysocks,
                                         libraries = "corncob")) %>%
    evaluate(list(n, wald, pbwald, lrt, pblrt))
  save(sim2, file = paste("TypeI_setting2_", nsim*100, "_", B, ".rda",
                          sep = ""))
} else if (arg1 == 3) {
  if (!file_test("-f", "draws3.rda")) {
    draws3 <- new_simulation(name = "setting3",
                             label = "Setting 3") %>%
      generate_model(make_betabin, seed = 1,
                     n = list(10,30,100), 
                     setting = 3,
                     vary_along = "n") %>%
      simulate_from_model(nsim = nsim, index = 1:100)
    save(draws3, file = "draws3.rda")
  }
  load("draws3.rda")
  sim3 <- draws3 %>%
    run_method(setting3, parallel = list(socket_names = mysocks,
                                         libraries = "corncob")) %>%
    evaluate(list(n, wald, pbwald, lrt, pblrt))
  save(sim3, file = paste("TypeI_setting3_", nsim*100, "_", B, ".rda",
                          sep = ""))
}



