abc_sim <- function (pars, data, tols, u, r) {
  tmp <- sim_SEIARD(pars, r)
  arr <- tmp[[1]]
  obs_sim <- tmp[[2]]
  if(any(is.na(arr[,,dim(arr)[3]]))) return(NA)
  
  met <- colSums(obs_sim) / N
  
  if (all(abs(met - data) <= tols)) {
    return(met)
  } else {
    return(NA)
  }
}

priors <- data.frame(
  parnames = c("alpha0", "alphar", "omega", "chi", "deltaI", "deltaA", "betaI", "eta", "pdet"),
  dist = rep("unif", 9),
  p1 = rep(0, 9),
  p2 = rep(0.5, 9),
  stringsAsFactors = F
)

u <- matrix(c(500000, 5, 0, 0, 0, 0,
              300000, 0, 0, 0, 0, 0,
              100000, 0, 0, 0, 0, 0), nrow=3, byrow = T)
data <- c(c1 = 0.0135378646213538, 
          c2 = 0.00520333333333333, 
          c3 = 0.00202)

abc_sim(pars, data, c(0.05, 0.05, 0.05), u, r)

post <- SimBIID::ABCSMC(
  x = data, 
  u = 1:3,
  priors = priors, 
  func = abc_sim, 
  # u = iniStates, 
  r = r,
  obs = obs,
  tols = c(c1 = 0.02, 
           c2 = 0.02, 
           c3 = 0.02), 
  ptol = 0.2, 
  ngen = 5, 
  npart = 200
)
