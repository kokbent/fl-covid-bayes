library(magrittr)
library(tidyverse)

states <- matrix(c(500000, 5, 0, 0, 0, 0,
                   300000, 0, 0, 0, 0, 0,
                   100000, 0, 0, 0, 0, 0), nrow=3, byrow = T)

N <- rowSums(states)
r <- matrix(c(0.9, 0.05, 0.01,
              0.05, 0.7, 0.01,
              0.01, 0.01, 0.5), nrow = 3, byrow = T)
Tt <- 300

params <- list()
params <- within(params,
                 {
                   alpha0 <- 0.0
                   alphar <- 0.35
                   omega <- 0.5
                   chi <- 0.3
                   deltaI <- 1/8
                   deltaA <- 1/6
                   betaI <- 0.7
                   betaA <- 1-betaI
                   eta <- 0.02
                   pdet <- 0.2
                 })

rmnom <- function (size, prob) sapply(size, function (x) rmultinom(1, x, prob))

calc_d <- function (states, N, params, r) {
  out <- with(
    params,
    {
      Alpha <- diag(alpha0, 3, 3) + alphar * r
      dS <- (- states[,1] / N * (Alpha %*% (states[,3] + omega * states[,4]))) %>%
        as.vector
      dE <- - dS - chi * states[,2]
      dI <- chi * betaI * states[,2] - (eta + deltaI) * states[,3]
      dA <- chi * betaA * states[,2] - deltaA * states[,4]
      dR <- deltaI * states[,3] + deltaA * states[,4]
      dD <- eta * states[,3]
      
      cbind(dS, dE, dI, dA, dR, dD)
    }
  )
  
  return(out)
}

calc_ds <- function (states, N, params, r) {
  out <- with(
    params,
    {
      cc <- length(N)
      Alpha <- diag(alpha0, 3, 3) + alphar * r
      dSo <- (states[,1] / N * (Alpha %*% (states[,3] + omega * states[,4]))) %>%
        as.vector
      dS <- - rpois(cc, dSo)
      
      dEo <- rbinom(cc, states[,2], chi)
      dE <- - dS - dEo
      
      dIi <- rbinom(cc, dEo, betaI)
      if(any(is.na(states[,3]))) {
        tmp <- matrix(NA, 3, 3)
      } else {
        tmp <- rmnom(states[,3], c(eta, deltaI, 1 - eta - deltaI))
      }
      
      dIoD <- tmp[1,]
      dIoR <- tmp[2,]
      dI <- dIi - dIoD - dIoR
      
      dAi <- dEo - dIi
      dAo <- rbinom(cc, states[,4], deltaA)
      dA <- dAi - dAo
      
      dR <- dIoR + dAo
      dD <- dIoD
      
      cbind(dS, dE, dI, dA, dR, dD, dIi)
    }
  )
  
  return(out)
}


update_states <- function(states, N, params, r) {
  d <- calc_ds(states, N, params, r)
  dIi <- d[,7]
  d <- d[,-7]
  states <- states + d
  
  return(list(states, dIi))
}

sim_SEIARD <- function (pars, r) {
  params <- within(params,
                   {
                     alpha0 <- pars[1]
                     alphar <- pars[2]
                     omega <- pars[3]
                     chi <- pars[4]
                     deltaI <- pars[5]
                     deltaA <- pars[6]
                     betaI <- pars[7]
                     eta <- pars[8]
                     pdet <- pars[9]
                   })
  
  states <- matrix(c(500000, 5, 0, 0, 0, 0,
                     300000, 0, 0, 0, 0, 0,
                     100000, 0, 0, 0, 0, 0), nrow=3, byrow = T)
  
  arr <- array(NA, c(3, 6, Tt+1),
               dimnames = list(NULL, c("S", "E", "I", "A", "R", "D")))
  arr[,,1] <- states
  obs <- matrix(NA, Tt+1, 3)
  obs[1,] <- rep(0, 3)
  for (i in 1:Tt) {
    if (i > 120) {
      tmp <- update_states(states, N, params, r/2)
    } else {
      tmp <- update_states(states, N, params, r)
    }
    # tmp <- update_states(states, N, params, r)
    dIi <- tmp[[2]]
    states <- tmp[[1]]
    obs[i+1,] <- rbinom(length(dIi), dIi, params$pdet)
    arr[,,i+1] <- states
  }

  return(list(arr, obs))
}

pars <- c(0, 0.35, 0.5, 0.3, 1/8, 1/6, 0.7, 0.02, 0.2)
mod <- sim_SEIARD(pars, r)
arr <- mod[[1]]
obs <- mod[[2]]

arr_df <- data.table::as.data.table(arr)
colnames(arr_df) <- c("county", "state", "time", "value")
arr_df$state <- factor(arr_df$state, levels = c("S", "E", "I", "A", "R", "D"))

ggplot(arr_df) +
  geom_line(aes(x=time, y=value, colour=as.factor(state))) +
  facet_wrap(~county, scales = "free_y")

obs_df <- data.table::as.data.table(obs)
colnames(obs_df) <- c("c1", "c2", "c3")
obs_df$time <- 1:nrow(obs_df)
obs_df <- pivot_longer(obs_df, -time, names_to = "county", values_to = "value")

ggplot(obs_df) +
  geom_line(aes(x=time, y=value/N, colour=county))