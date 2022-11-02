library(tidyverse)
library(dplyr)
library(ggplot2)
library(pROC)

# set up -----------------------------------------------------------------------
set.seed(1)
n <- 1000

# generate risk_t1
risk_t1 <-  runif(n, 0, 1)

# intervention when risk_t1 > cmi_cutoff
cmi_cutoff <- 0.5

# intervention effectiveness/probability of the intervention works (fixed, varies in different trials)
b_intervention <- 1 # set to 1  for null case

# parameter to find the appropriate AUC/AUC = 0.80
beta_1 <- 0.9

# generate unintervened outcome from risk_t1 probabilistically------------------
gen_outcome1 <- function(risk_t1, n = n, beta_1){
  outcome1 <- rbinom(n, 1, beta_1 * risk_t1) # beta_1 vary to vary AUC
  return (as.data.frame(cbind(risk_t1, outcome1)))
}
df
a <- gen_outcome1(risk_t1, n = n, beta_1)
auc(a[,'outcome1'], a[,'risk_t1'])

# generate cmi (binary) depending on cmi_cutoff---------------------------------
gen_cmi <- function(risk_t1, n = n){
  cmi <- case_when(risk_t1 > cmi_cutoff ~ 1, TRUE ~ 0)
  return (as.matrix(cmi))
}

# generate risk_t2 depending on cmi & risk_t1-----------------------------------
gen_risk_t2 <- function(risk_t1, n = n, b_intervention, beta_1){
  cmi <- gen_cmi(risk_t1, n)
  outcome_t1 <- gen_outcome1(risk_t1, n, beta_1)
  r <- cbind(outcome_t1, cmi)

  # if cmi == 1 & outcome_no_interventnion == 1 then risk_t2 = risk_t1 * b_intervention
  risk_t2 <- case_when((r[,3] == 1 & r[,2] == 1) ~ r[,1] * b_intervention,
                       TRUE ~ risk_t1)

  risks <- cbind(r, risk_t2)
  colnames(risks) <- c("risk_t1", "outcome_no_intervention", "intervention", "risk_t2")
  return (as.matrix(risks))
}

b <- gen_risk_t2(risk_t1, n = n, b_intervention, beta_1)
# table(b[,2], b[,3])

# generate outcome from risk_t2-------------------------------------------------
gen_outcome <- function(risk_t1, n = n, b_intervention, beta_1){
  risks <- as.data.frame(gen_risk_t2(risk_t1, n, b_intervention, beta_1))

  outcome_no_intervention <- risks[,2]
  risk_t2 <- risks[,4]

  # generate all observed outcome probabilistically
  # oo: if outcome_no_intervention = 0 & intervention = 0 --> old outcome
  # ol: if outcome_no_intervention = 0 & intervention = 1 --> old outcome
  # lo: if outcome_no_intervention = 1 & intervention = 0 --> old outcome
  # ll: if outcome_no_intervention = 1 & intervention = 1 --> rbinom(n, 1, risks_t2)
  outcome_intervention <- rbinom(n, 1, risk_t2)
  oo <- which((risks[,2] == 0 & risks[,3] == 0))
  ol <- which(risks[,2] == 0 & risks[,3] == 1)
  lo <- which(risks[,2] == 1 & risks[,3] == 0)
  outcome_intervention[oo] <- outcome_no_intervention[oo]
  outcome_intervention[ol] <- outcome_no_intervention[ol]
  outcome_intervention[lo] <- outcome_no_intervention[lo]

  return (as.matrix(cbind(risks, outcome_intervention)))
}


# model training-----------------------------------------------------------------
outcomes <- as.data.frame(gen_outcome(risk_t1, n = n, b_intervention, beta_1))

mean(outcomes$risk_t1)
mean(outcomes$risk_t2)
mean(outcomes$intervention)
mean(outcomes$outcome_no_intervention)
mean(outcomes$outcome_intervention)

x <- outcomes[,c(1,3)]
y <-  outcomes[,5]

# train/test split
train <- sample(1:nrow(x), nrow(x)/2)
x.train <- x[train,]
y.train <- y[train]

test <- (-train)
x.test <- x[test,]
y.test <- y[test]

# fit logistic regression
# Difference method

# indirect effect
fit <- glm(y.train ~ ., data = x.train, family = "binomial")
# coefficient of risk_t1: 6.8678
# indirect effect = 6.8678-4.4338 = 2.434

# direct effect
fit_2 <- glm(y.train ~ ., data = as.data.frame(x.train[,1]), family = "binomial")
# coefficient of risk_t1: 4.4338

# Evaluate-----------------------------------------------------------------------
y.pred <- predict(fit, x.test, type = "response")

# AUC
AUC <- auc(y.test, y.pred)
AUC

# AUC grid----------------------------------------------------------------------
AUC.grid <- rep(NA, 10)
for (i in seq(from=0.1, to=1, by=0.1)) {

  outcomes <- as.data.frame(gen_outcome(risk_t1,
                                        n,
                                        b_intervention = i,
                                        beta_1))
  x <- outcomes[,c(1,3)]
  y <-  outcomes[,5]
  train <- sample(1:nrow(x), nrow(x)/2)
  x.train <- x[train,]
  y.train <- y[train]
  test <- (-train)
  x.test <- x[test,]
  y.test <- y[test]
  fit <- glm(y.train ~ ., data = x.train, family = "binomial")
  y.pred <- predict(fit, x.test, type = "response")
  AUC.grid[i*10] <- auc(y.test, y.pred)
}

plot(AUC.grid, xlab = "b_intervention * 10", ylab = "AUC", type = "l")


# Notes-------------------------------------------------------------------------
# % outcome changes from 1 to 0 after intervention?
changed <- which(outcomes[,2] == 1 & outcomes[,5] == 0)
unchanged <- which(outcomes[,2] == 1 & outcomes[,5] == 1)
length(changed)/(length(unchanged) + length(changed))


# When the effectiveness of intervention is stronger we expect AUC to decrease.
# At lower b_intervention,the model might predict everything to be 0 causing a high
# AUC?

# next week--------------------------------------------------------------------
# vary b_intervention (many simulation)
# vary percentage of cmi/cmi cutoff



