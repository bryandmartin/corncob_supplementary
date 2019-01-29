## @knitr models

make_betabin <- function(n, cons, setting) {
  # n is sample size
  # Setting is 1 or 2
  
  # Helpful: paste0(collapse = ",")
  fullM <- c(57183,51769,24557,50768,10990,35795,21280,20117,53869,33660,18457,
             24409,23141,58655,16394,43292,22485,14635,26561,21158,18429,19058,
             17504,20540,13018,25244,17531,16409,20010,16007,7821)
  
  X <- X_star <- matrix(seq(0,1,length = n), nrow = n, ncol = 1)
  colnames(X) <- paste0("X", 1)
  colnames(X_star) <- paste0("Xst", 1)
  
  beta <- c(-5.16655052974227,-2.46370460866063)
  beta_star <- c(-5.13423469699394,-3.87644320689491)
  
  if (setting == 1) {
    beta[2] <- beta[2]*cons
  } else if (setting == 2) {
    beta_star[2] <- beta_star[2]*cons
  }
  
  # Generate mu_i where logit(mu) = X*beta
  mu_i <- invlogit(cbind(1, X) %*% beta)
  
  # Generate phi_i where fishZ(phi) = X_star * beta_star
  phi_i <- invlogit(cbind(1, X_star) %*% beta_star)
  
  # Draw n with replacement
  M_i <- sample(fullM, n, replace = TRUE)
  

  new_model(name = "betabinomial-model",
            label = sprintf("Betabinomial (n = %s, setting = %s, cons = %s)",
                            n, setting, cons),
            params = list(n = n,
                          M = M_i,
                          mu = mu_i,
                          phi = phi_i,
                          X = X,
                          Xstar = X_star,
                          cons = cons),
            simulate = function(n, M, mu, phi, nsim) {
              W_i <- replicate(nsim, 
                               VGAM::rbetabinom(n = n, size = M, prob = mu,
                                                rho = phi), simplify = FALSE)
            })
}


