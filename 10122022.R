library(tidyverse)
library(dplyr)
library(ggplot2)
library(pROC)

# set up ------------------------------------------------------------------
set.seed(1)
n <- 1000

# generate risk_t1
risk_t1 <-  runif(n, 0, 1)


# threshold of detecting outcome
threshold <- 0.6

# intervention effectiveness (fixed)
b_intervention <- 0.75


# generate cmi (binary)
# mechanism 1: probabilistic
gen_cmi <- function(risk_t1, beta_rt, n = n){
  l <- beta_rt * risk_t1
  pr <- 1/(1 + exp(-l))
  cmi <- rbinom(n, 1, pr)
  return (as.matrix(cbind(risk_t1, pr, cmi)))
}

# different mechanism to receive intervention (i.e. intervene if risk_t1 > some threshold)

# generate risk_t2
gen_risk_t2 <- function(risk_t1, beta_rt, n = n, b_intervention){
  cmi <- gen_cmi(risk_t1, beta_rt, n)
  risks <- cbind(cmi, cmi[,1])

  # if intervention is involved, vary the effectiveness of the intervention
  # percent risk left after intervention
  intervention <- which(risks[,3] == 1)

  # risk_t2: risks[,4]
  # relative risk: exponential
  risks[,4][intervention] <- risks[,1] * b_intervention # effectiveness (vary) of intervention; treatment effect should be fixed; relative risk

  # old problem persist? 0.631420228, 1, 0.738821418
  return (as.matrix(cbind(risks)))
}

# generate outcome
gen_outcome <- function(risk_t1, beta_rt, n = n, b_intervention){
  risks <- gen_risk_t2(risk_t1, beta_rt, n = n, b_intervention)
  risks <- cbind(risks, 0)
  risk_t2 <- risks[,4]

  # outcome: risks[,6]
  risks[,5] <- rbinom(n, 1, risk_t2)
}


# how strong intervention needs to be for outcome to change? (See pic)