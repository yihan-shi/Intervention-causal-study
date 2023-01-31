# load libraries ----------------------------------------------------------
library(tidyverse)
library(dplyr)
library(ggplot2)
library(pROC)
library(logbin)
library(epitools)

# generate data -----------------------------------------------------------
set.seed(1)
n <- 20000
repeats <- 20 # for simulation
cmi_cutoff <- 0.7
b_intervention <- 0.2

risk_t1 <- runif(n, 0, 1)
prevalence <- c(-6, -4, -2.5)
prev <- c("low", "mid", "high")

gen_outcome1 <- function(risk_t1, n = n, prevalence){
  xb <- prevalence + 5 * risk_t1
  p <- exp(xb)/(1 + exp(xb))
  outcome_noint <- rbinom(n, 1, p)
  df <- as.data.frame(cbind(risk_t1, outcome_noint))
  return (df)
}

a <- gen_outcome1(risk_t1, n, prevalence = -2.5)

gen_cmi <- function(risk_t1, n = n, cmi_cutoff, prevalence){
  df <- gen_outcome1(risk_t1, n, prevalence)
  df <- df %>%
    mutate(cmi = ifelse(risk_t1 > cmi_cutoff, 1, 0)) # the higher the cutoff, the less people getting intervention

  return (df)
}

gen_risk_t2 <- function(risk_t1, n = n, b_intervention, cmi_cutoff, prevalence){
  # if cmi == 1 & outcome_noint == 1, then risk_t2 = risk_t1 * b_intervention
  df <- gen_cmi(risk_t1, n, cmi_cutoff, prevalence)
  df <- df %>%
    mutate(risk_t2 = ifelse((outcome_noint == 1 & cmi == 1),
                            risk_t1 * b_intervention,
                            risk_t1))
  return (df)
}

gen_outcome_f1 <- function(risk_t1, n = n, b_intervention,cmi_cutoff, prevalence){
  df <- gen_risk_t2(risk_t1, n, b_intervention, cmi_cutoff, prevalence)
  condition <- df %>%
    filter(outcome_noint == 1 & cmi == 1)

  # generate all observed outcome probabilistically
  # oo: if outcome_no_intervention = 0 & intervention = 0 --> old outcome
  # ol: if outcome_no_intervention = 0 & intervention = 1 --> old outcome
  # lo: if outcome_no_intervention = 1 & intervention = 0 --> old outcome
  # ll: if outcome_no_intervention = 1 & intervention = 1 --> rbinom(n, 1, risks_t2)

  num_obs <- dim(condition)[1]

  xb <- prevalence + 5 * df$risk_t2
  p <- exp(xb)/(1 + exp(xb))

  df_final <- df %>%
    mutate(outcome_int = ifelse((outcome_noint  == 1 & cmi == 1),
                                # rbinom(num_obs, 1, df$risk_t2),
                                rbinom(num_obs, 1, p),
                                outcome_noint))

  return (df_final)
}

# load dataset ------------------------------------------------------------
df <- as.data.frame(gen_outcome_f1(risk_t1,
                                  n,
                                  b_intervention = 0.1,
                                  cmi_cutoff = 0.5,
                                  prevalence = -4))

x <- df %>%
  select(-outcome_int)
y <- df %>%
  select(outcome_int)

# train/test split
#train <- sample(1:nrow(x), nrow(x)/1.2)
train <- 1:round(0.8 * nrow(x))
x.train <- x[train,]
y.train <- y[train,]

#test <- (-train)
test <- (round(0.8 * nrow(x)) + 1):nrow(x)
x.test <- x[test,]
y.test <- y[test,]

# Intervention type 1: Resource-constrained (cutoff approach)

# Correction 1: ---------------------------------------------------------------
# Use log(p) = intercept + \beta_{i} * X_i + \beta_{r} * X_r to estimate \beta_{i}, which is the
# effective of intervention. New risk = Initial risk * \beta_{i}

fit_lb <- logbin(y.train ~ risk_t1, data = as.data.frame(x.train), method = "em") # method determines which algorithm to use to find the MLE
car::S(fit_lb)

y.pred_resp_lb <- predict(fit_lb, as.data.frame(x.test), type = "response")
auc(y.test, y.pred_resp)

# known: risk_t1, CMI
fit_lb2 <- logbin(y.train ~ risk_t1 + cmi, data = as.data.frame(x.train), method = "em", maxit = 50000) # method determines which algorithm to use to find the MLE
car::S(fit_lb2)
y.pred_resp_lb <- predict(fit_lb2, as.data.frame(x.test), type = "response")
auc(y.test, y.pred_resp_lb)


# Correction 2:  ---------------------------------------------------------------
# Stratified analysis: calculate weighted AUC for CMI = 1 and CMI = 0.
# Finding:
# Weighted AUC improves

set.seed(10)
model_sum <- expand.grid(prev = prevalence,
                         b_int = c(0.1, 0.5, 0.9),
                         cmi = c(0.1, 0.5, 0.9),
                         mean_squared_error = NA)
row <- 1
for (row in 1:dim(model_sum)[1]){
  # mse <- rep(NA, 100)
  auc <- 0
  auc_uncorrected <- 0

  for (i in 1:100){
    df <- as.data.frame(gen_outcome_f1(risk_t1,
                                      n,
                                      b_intervention = model_sum$b_int[row],
                                      cmi_cutoff = model_sum$cmi[row],
                                      prevalence = model_sum$prev[row]))

    fit <- glm(outcome_int ~ risk_t1, data = df, family = binomial)
    y.pred <- predict(fit, as.data.frame(df), type = "response")

    have_cmi <-subset(df, cmi == 1)
    no_cmi <-subset(df, cmi == 0)

    fit_cmi <- glm(outcome_int ~ risk_t1, data = have_cmi, family = binomial)
    y.pred_cmi <- predict(fit_cmi, as.data.frame(have_cmi), type = "response")
    fit_no_cmi <- glm(outcome_int ~ risk_t1, data = no_cmi, family = binomial)
    y.pred_no_cmi <- predict(fit_no_cmi, as.data.frame(no_cmi), type = "response")

    auc = auc + dim(have_cmi)[1]/dim(df)[1] * auc(have_cmi$outcome_int, y.pred_cmi) + dim(no_cmi)[1]/dim(df)[1] * auc(no_cmi$outcome_int, y.pred_no_cmi)
    auc_uncorrected = auc_uncorrected + auc(df$outcome_int, y.pred)
  }

  model_sum$auc[row] <- auc/100
  model_sum$auc_uncorrected[row] <- auc_uncorrected/100
  print(paste0("row", row))
}

model_sum <- model_sum %>%
  arrange(prev) %>%
  mutate(auc_diff = auc - auc_uncorrected)

# ---------------------------------------------------------------------------
have_cmi <-subset(df, cmi == 1)

x <- have_cmi %>%
  select(-outcome_int)
y <- have_cmi %>%
  select(outcome_int)

train <- 1:round(0.7 * nrow(x))
x.train <- x[train,]
y.train <- y[train,]

# test <- (-train)
test <- (round(0.7 * nrow(x)) + 1):nrow(x)
x.test <- x[test,]
y.test <- y[test,]

x.test <- x.test %>%
  filter(risk_t1 > range(x.train$risk_t1)[1] | risk_t1 < range(x.train$risk_t1)[2])

fit_havecmi <- logbin(y.train ~ risk_t1, data = as.data.frame(x.train), method = "em")
car::S(fit_havecmi)
# y.pred <- predict(fit_havecmi, as.data.frame(x.test), type = "response")
# auc(y.test, y.pred_resp)

# ---------------------------------------------------------------------------
no_cmi <-subset(df, cmi == 0)

x <- no_cmi %>%
  select(-outcome_int)
y <- no_cmi %>%
  select(outcome_int)

train <- 1:round(0.7 * nrow(x))
x.train <- x[train,]
y.train <- y[train,]

# test <- (-train)
test <- (round(0.7 * nrow(x)) + 1):nrow(x)
x.test <- x[test,]
y.test <- y[test,]

fit_nocmi <- logbin(y.train ~ risk_t1, data = as.data.frame(x.train), method = "em") # method determines which algorithm to use to find the MLE
car::S(fit_nocmi)
# y.pred <- predict(fit_nocmi, as.data.frame(x.test), type = "response")
# auc(y.test, y.pred_resp)


# # Intervention type 2: Probabilistic approach

# resource ----------------------------------------------------------------
# https://www.statulator.com/blog/conducting-stratified-analyses/
# https://bookdown.org/rwnahhas/RMPH/blr-log-binomial.html
# https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_multivariable/bs704_multivariable3.html#:~:text=The%20Cochran%2DMantel%2DHaenszel%20method%20is%20a%20technique%20that%20generates,and%20a%20dichotomous%20risk%20factor.

