## @knitr models

make_betabin <- function(n, setting, beta, beta_star) {
  # n is sample size
  # Setting is 1, 2, or 3
  
  # Helpful: paste0(collapse = ",")
  fullM <- c(57183,51769,24557,50768,10990,35795,21280,20117,53869,33660,18457,
             24409,23141,58655,16394,43292,22485,14635,26561,21158,18429,19058,
             17504,20540,13018,25244,17531,16409,20010,16007,7821)
  
  X <- X_star <- matrix(seq(0,1,length = n), nrow = n, ncol = 1)
  colnames(X) <- paste0("X", 1)
  colnames(X_star) <- paste0("Xst", 1)
  
  if (setting == 1) {
    beta <- c(-5.7537429091197,0)
    beta_star <- c(-5.24001799394086,0)
  } else if (setting == 2) {
    beta <- c(-5.3578960694989,-1.11962916625361)
    beta_star <- c(-5.69049638516877,0)
  } else {
    beta <- c(-5.51303444999736,0)
    beta_star <- c(-5.38197099149681,0.700660697997067)
  }
  
  # Generate mu_i where logit(mu) = X*beta
  mu_i <- invlogit(cbind(1, X) %*% beta)
  
  # Generate phi_i where fishZ(phi) = X_star * beta_star
  phi_i <- invlogit(cbind(1, X_star) %*% beta_star)
  
  # Draw n with replacement
  M_i <- sample(fullM, n, replace = TRUE)
  

  new_model(name = "betabinomial-model",
            label = sprintf("Betabinomial (n = %s, setting = %s)", n, setting),
            params = list(n = n,
                          M = M_i,
                          mu = mu_i,
                          phi = phi_i,
                          X = X,
                          Xstar = X_star),
            simulate = function(n, M, mu, phi, nsim) {
              W_i <- replicate(nsim, 
                               VGAM::rbetabinom(n = n, size = M, prob = mu,
                                                rho = phi), simplify = FALSE)
            })
}


