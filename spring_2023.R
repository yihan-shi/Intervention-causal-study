# load libraries ----------------------------------------------------------
library(tidyverse)
library(dplyr)
library(ggplot2)
library(pROC)
library(logbin)
library(epitools)

# generate data -----------------------------------------------------------
set.seed(10)
n <- 40000
repeats <- 50 # for simulation
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


# Intervention type 1: Resource-constrained (cutoff approach)
# Correction 1: ---------------------------------------------------------------
# Use log(p) = intercept + \beta_{i} * X_i + \beta_{r} * X_r to estimate \beta_{i}, which is the
# effective of intervention. New risk = Initial risk * \beta_{i}
# Finding:
# Most weighted AUC improve (some errors)

model_sum_lb <- expand.grid(prev = prevalence,
                         b_int = seq(0.1, 0.9, 0.1),
                         cmi = seq(0.1, 0.9, 0.1),
                         auc = NA,
                         auc_uncorrected = NA)
                         # b_int = c(0.1, 0.5, 0.9),
                         # cmi = c(0.1, 0.5, 0.9))

set.seed(10)
for (row in 1:dim(model_sum_lb)[1]){
# for (row in 111:111) {
  print(paste0("row", row))
  auc <- 0
  auc_uncorrected <- 0

  for (i in 1:repeats){
  # for (i in 1:6){
    # load data
    df <- as.data.frame(gen_outcome_f1(risk_t1,
                                       n,
                                       b_intervention = model_sum_lb$b_int[row],
                                       cmi_cutoff = model_sum_lb$cmi[row],
                                       prevalence = model_sum_lb$prev[row]))

    x <- df %>%
      select(-c(outcome_int, risk_t2, outcome_noint))
    y <- df %>%
      select(outcome_int)

    # train/test split
    train <- 1:round(0.8 * nrow(x))
    x.train <- x[train,]
    y.train <- y[train,]
    test <- (round(0.8 * nrow(x)) + 1):nrow(x)
    x.test <- x[test,]
    y.test <- y[test,]

    # original auc
    # ERROR HANDLING
    fit <- glm(y.train ~ risk_t1, data = x.train, family = binomial(link = "logit"))

    possibleError <- tryCatch(
      glm(y.train ~ risk_t1, data = x.train, family = binomial(link = "log"), start = coef(fit)),
      error=function(e) e
      )

    if(inherits(possibleError, "error")) {
      auc = 0
      print("supply starting value error")
      break
    }
    fit_again <- glm(y.train ~ risk_t1, data = x.train, family = binomial(link = "log"), start = coef(fit))

    # tapply(x.train$risk_t1, y.train, summary)
    # summary(fit)

    x.test_1 <- x.test %>%
      select(risk_t1)
    y.pred_resp <- predict.glm(fit_again, x.test_1, type = "response")
    auc_uncorrected = auc_uncorrected + auc(y.test, y.pred_resp)

    # fit model with new information: CMI
    # ERROR HANDLING
    fit_lb <- glm(y.train ~ risk_t1 + as.factor(cmi), data = x.train, family = binomial(link = "logit"))
    possibleError <- tryCatch(
      glm(y.train ~ risk_t1 + as.factor(cmi), data = x.train, family = binomial(link = "log"), start = coef(fit_lb)),
      error=function(e) e
    )

    if(inherits(possibleError, "error")) {
      print("log-binomial model has invalid starting values")
      auc = 0
      break
    }
    fit_lb_again <- glm(y.train ~ risk_t1 + as.factor(cmi), data = x.train, family = binomial(link = "log"), start = coef(fit_lb))

    # deterministically can't compare cmi at fixed level of risk_t1, problematic method bc need to supply start value

    # for the patients with observation, use exp(\beta_{i}) to adjust their risk
    x.train_2 <- x.train %>%
      mutate(new_risk = ifelse(cmi == 1, risk_t1 * exp(coef(fit_lb_again)[3]),
                               risk_t1))
    x.test_2 <- x.test %>%
      mutate(new_risk = ifelse(cmi == 1, risk_t1 * exp(coef(fit_lb_again)[3]),
                               risk_t1))

    # ERROR HANDLING
    fit_lb2 <- glm(y.train ~ new_risk, data = x.train_2, family = binomial(link = "logit"))
    possibleError <- tryCatch(
      glm(y.train ~ new_risk, data = x.train_2, family = binomial(link = "log"), start = coef(fit_lb2)),
      error=function(e) e
    )

    if(inherits(possibleError, "error")) {
      print("log-binomial model has invalid starting values")
      auc = 0
      break
    }
    fit_lb2_again <- glm(y.train ~ new_risk, data = x.train_2, family = binomial(link = "log"), start = coef(fit_lb2))


    # use new_risk to fit model
    x.test_2 <- x.test_2 %>%
      select(new_risk)
    y.pred_resp_2 <- predict.glm(fit_lb2, x.test_2, type = "response")
    auc = auc + auc(y.test, y.pred_resp_2)
  }
  model_sum_lb$auc[row] <- auc/repeats
  model_sum_lb$auc_uncorrected[row] <- auc_uncorrected/repeats
}

model_sum_lb <- model_sum_lb %>%
  arrange(prev) %>%
  mutate(auc_diff = auc - auc_uncorrected)

write.csv(model_sum_lb, "model_sum_lb_0216.csv")

# Correction 2:  ---------------------------------------------------------------
# Stratified analysis: calculate weighted AUC for CMI = 1 and CMI = 0.
# Finding:
# Weighted AUC didn't improve consistently (many even decreased?)

set.seed(10)
model_sum <- expand.grid(prev = prevalence,
                            b_int = seq(0.1, 0.9, 0.1),
                            cmi = seq(0.1, 0.9, 0.1),
                            auc = NA,
                            auc_uncorrected = NA)

# row <- 1
for (row in 1:dim(model_sum)[1]){
# for (row in 1:1){
  print(paste0("row", row))
  # mse <- rep(NA, 100)
  auc <- 0
  auc_uncorrected <- 0

  for (i in 1:repeats){
  # for (i in 20:20) {
    df <- as.data.frame(gen_outcome_f1(risk_t1,
                                      n,
                                      b_intervention = model_sum$b_int[row],
                                      cmi_cutoff = model_sum$cmi[row],
                                      prevalence = model_sum$prev[row]))

    fit <- glm(outcome_int ~ risk_t1, data = df, family = binomial)
    y.pred <- predict(fit, as.data.frame(df), type = "response")

    have_cmi <- df %>%
      filter(cmi == 1)
    no_cmi <- df %>%
      filter(cmi == 0)

    fit_cmi <- glm(outcome_int ~ risk_t1, data = have_cmi, family = binomial)
    y.pred_cmi <- predict(fit_cmi, as.data.frame(have_cmi), type = "response")
    fit_no_cmi <- glm(outcome_int ~ risk_t1, data = no_cmi, family = binomial)
    y.pred_no_cmi <- predict(fit_no_cmi, as.data.frame(no_cmi), type = "response")

    possibleError <- tryCatch(
      auc(no_cmi$outcome_int, y.pred_no_cmi),
      auc(have_cmi$outcome_int, y.pred_cmi),
      auc(df$outcome_int, y.pred),
      error=function(e) e
    )

    if(inherits(possibleError, "error")) {
      auc = 0
      print("No-intervention group have no outcome")
      break
    }

    auc = auc + dim(have_cmi)[1]/dim(df)[1] * auc(have_cmi$outcome_int, y.pred_cmi) + dim(no_cmi)[1]/dim(df)[1] * auc(no_cmi$outcome_int, y.pred_no_cmi)
    auc_uncorrected = auc_uncorrected + auc(df$outcome_int, y.pred)
  }

  model_sum$auc[row] <- auc/repeats
  model_sum$auc_uncorrected[row] <- auc_uncorrected/repeats
}

model_sum <- model_sum %>%
  arrange(prev) %>%
  mutate(auc_diff = auc - auc_uncorrected)

write.csv(model_sum, "model_sum_0216.csv")

# Intervention type 2: Probabilistic approach ---------------------------------
# new functions
gen_cmi_2 <- function(risk_t1, n = n, prevalence){
  df <- gen_outcome1(risk_t1, n, prevalence)

  xb <- prevalence + 5 * risk_t1
  p <- exp(xb)/(1 + exp(xb))

  df <- df %>%
    mutate(cmi = rbinom(n, 1, p),
           prob = p) # risk is proportional to the likelihood of cmi;adjust for auc & fraction of ppl getting intervention
  return (df)
}

df <- gen_cmi_2(risk_t1, n, prevalence = -2.5)
a <- df %>%
  filter(cmi == 1)
mean(a$prob)
mean(a$risk)

b <- df %>%
  filter(cmi == 0)
mean(b$risk)

gen_risk_t2_2 <- function(risk_t1, n = n, b_intervention, prevalence){
  # if cmi == 1 & outcome_noint == 1, then risk_t2 = risk_t1 * b_intervention
  df <- gen_cmi_2(risk_t1, n, prevalence)
  df <- df %>%
    mutate(risk_t2 = ifelse((outcome_noint == 1 & cmi == 1),
                            risk_t1 * b_intervention,
                            risk_t1))
  return (df)
}

gen_outcome_f2 <- function(risk_t1, n = n, b_intervention,prevalence){
  df <- gen_risk_t2_2(risk_t1, n, b_intervention, prevalence)
  condition <- df %>%
    filter(outcome_noint == 1 & cmi == 1)

  num_obs <- dim(condition)[1]

  xb <- prevalence + 5 * df$risk_t2
  p <- exp(xb)/(1 + exp(xb))

  df_final <- df %>%
    mutate(outcome_int = ifelse((outcome_noint  == 1 & cmi == 1),
                                rbinom(num_obs, 1, p),
                                outcome_noint))

  return (df_final)
}

df <- as.data.frame(gen_outcome_f2(risk_t1,
                                   n,
                                   b_intervention = 0.1,
                                   prevalence = -2.5))
mean(df$cmi)
mean(df$outcome_noint)
mean(df$outcome_int)

# Intervention type 2: Probabilistic approach ---------------------------------
# Correction 1: ---------------------------------------------------------------
# Use log(p) = intercept + \beta_{i} * X_i + \beta_{r} * X_r to estimate \beta_{i}, which is the
# effective of intervention. New risk = Initial risk * \beta_{i}
# problem: many "supply starting value error"


model_sum_lb2 <- expand.grid(prev = prevalence,
                            b_int = seq(0.1, 0.9, 0.1),
                            auc = NA,
                            auc_uncorrected = NA)

set.seed(10)
for (row in 1:dim(model_sum_lb2)[1]){
# for (row in 1:1){
  print(paste0("row", row))
  auc <- 0
  auc_uncorrected <- 0

  for (i in 1:repeats){
  # for (i in 1:1){
    # load data
    df <- as.data.frame(gen_outcome_f2(risk_t1,
                                       n,
                                       b_intervention = model_sum_lb2$b_int[row],
                                       prevalence = model_sum_lb2$prev[row]))

    x <- df %>%
      select(-c(outcome_int, risk_t2, outcome_noint))
    y <- df %>%
      select(outcome_int)

    # train/test split
    train <- 1:round(0.8 * nrow(x))
    x.train <- x[train,]
    y.train <- y[train,]
    test <- (round(0.8 * nrow(x)) + 1):nrow(x)
    x.test <- x[test,]
    y.test <- y[test,]

    # original auc
    fit <- glm(y.train ~ risk_t1, data = x.train, family = binomial(link = "logit"))
    summary(fit)
    # ERROR HANDLING
    possibleError <- tryCatch(
      glm(y.train ~ risk_t1, data = x.train, family = binomial(link = "log"), start = coef(fit)),
      error=function(e) e
    )

    if(inherits(possibleError, "error")) {
      auc = 0
      print("supply starting value error")
      break
    }
    fit_again <- glm(y.train ~ risk_t1, data = x.train, family = binomial(link = "log"), start = coef(fit))

    x.test_1 <- x.test %>%
      select(risk_t1)
    y.pred_resp <- predict.glm(fit_again, x.test_1, type = "response")
    auc_uncorrected = auc_uncorrected + auc(y.test, y.pred_resp)

    # fit model with new information: CMI
    fit_lb <- glm(y.train ~ risk_t1 + as.factor(cmi), data = x.train, family = binomial(link = "logit"))
    # ERROR HANDLING
    possibleError <- tryCatch(
      glm(y.train ~ risk_t1 + as.factor(cmi), data = x.train, family = binomial(link = "log"), start = coef(fit_lb)),
      error=function(e) e
    )

    if(inherits(possibleError, "error")) {
      auc = 0
      break
    }
    fit_lb_again <- glm(y.train ~ risk_t1 + as.factor(cmi), data = x.train, family = binomial(link = "log"), start = coef(fit_lb))

    # deterministically can't compare cmi at fixed level of risk_t1, problematic method bc need to supply start value

    # for the patients with observation, use exp(\beta_{i}) to adjust their risk
    x.train_2 <- x.train %>%
      mutate(new_risk = ifelse(cmi == 1, risk_t1 * exp(coef(fit_lb_again)[3]),
                               risk_t1))
    x.test_2 <- x.test %>%
      mutate(new_risk = ifelse(cmi == 1, risk_t1 * exp(coef(fit_lb_again)[3]),
                               risk_t1))


    fit_lb2 <- glm(y.train ~ new_risk, data = x.train_2, family = binomial(link = "logit"))
    # ERROR HANDLING
    possibleError <- tryCatch(
      glm(y.train ~ new_risk, data = x.train_2, family = binomial(link = "log"), start = coef(fit_lb2)),
      error=function(e) e
    )

    if(inherits(possibleError, "error")) {
      auc = 0
      print("supply starting value error")
      break
    }
    fit_lb2_again <- glm(y.train ~ new_risk, data = x.train_2, family = binomial(link = "log"), start = coef(fit_lb2))

    # use new_risk to fit model
    x.test_2 <- x.test_2 %>%
      select(new_risk)
    y.pred_resp_2 <- predict.glm(fit_lb2, x.test_2, type = "response")
    auc = auc + auc(y.test, y.pred_resp_2)
  }
  model_sum_lb2$auc[row] <- auc/repeats
  model_sum_lb2$auc_uncorrected[row] <- auc_uncorrected/repeats
}

model_sum_lb2 <- model_sum_lb2 %>%
  arrange(prev) %>%
  mutate(auc_diff = auc - auc_uncorrected)

write.csv(model_sum_lb2, "model_sum_lb2_0216.csv")

# Correction 2:  ---------------------------------------------------------------
# Stratified analysis: calculate weighted AUC for CMI = 1 and CMI = 0.
# Finding:
# Weighted AUC improves

model_sum2 <- expand.grid(prev = prevalence,
                         b_int = seq(0.1, 0.9, 0.1),
                         auc = NA,
                         auc_uncorrected = NA)

set.seed(10)
for (row in 1:dim(model_sum2)[1]){
# for (row in 243:243){
  print(paste0("row", row))
  auc <- 0
  auc_uncorrected <- 0

  for (i in 1:repeats){
    df <- as.data.frame(gen_outcome_f2(risk_t1,
                                       n,
                                       b_intervention = model_sum2$b_int[row],
                                       prevalence = model_sum2$prev[row]))

    fit <- glm(outcome_int ~ risk_t1, data = df, family = binomial)
    y.pred <- predict(fit, as.data.frame(df), type = "response")

    have_cmi <- df %>%
      filter(cmi == 1)
    no_cmi <- df %>%
      filter(cmi == 0)

    fit_cmi <- glm(outcome_int ~ risk_t1, data = have_cmi, family = binomial)
    y.pred_cmi <- predict(fit_cmi, as.data.frame(have_cmi), type = "response")
    fit_no_cmi <- glm(outcome_int ~ risk_t1, data = no_cmi, family = binomial)
    y.pred_no_cmi <- predict(fit_no_cmi, as.data.frame(no_cmi), type = "response")

    possibleError <- tryCatch(
      auc(no_cmi$outcome_int, y.pred_no_cmi),
      auc(have_cmi$outcome_int, y.pred_cmi),
      auc(df$outcome_int, y.pred),
      error=function(e) e
    )

    if(inherits(possibleError, "error")) {
      auc = 0
      print("No-intervention group have no outcome")
      break
    }

    auc = auc + dim(have_cmi)[1]/dim(df)[1] * auc(have_cmi$outcome_int, y.pred_cmi) + dim(no_cmi)[1]/dim(df)[1] * auc(no_cmi$outcome_int, y.pred_no_cmi)
    auc_uncorrected = auc_uncorrected + auc(df$outcome_int, y.pred)
  }

  model_sum2$auc[row] <- auc/repeats
  model_sum2$auc_uncorrected[row] <- auc_uncorrected/repeats
}

model_sum2 <- model_sum2 %>%
  arrange(prev) %>%
  mutate(auc_diff = auc - auc_uncorrected)

write.csv(model_sum2, "model_sum2_0216.csv")


# resource ----------------------------------------------------------------
# https://www.statulator.com/blog/conducting-stratified-analyses/
# https://bookdown.org/rwnahhas/RMPH/blr-log-binomial.html
# https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_multivariable/bs704_multivariable3.html#:~:text=The%20Cochran%2DMantel%2DHaenszel%20method%20is%20a%20technique%20that%20generates,and%20a%20dichotomous%20risk%20factor.


# ---
# true_auc, corrected_1, corrected_2, uncorrected (side-by-side boxplot/violin plot)
# show auc plot under 2 kinds of intervention

# process results --------------------------------------------------------------
# intervention 1
i1 <- read.csv("model_sum_0216.csv")
i1_lb <- read.csv("model_sum_lb_0216.csv")

i1 <- i1 %>%
  filter(auc != 0) %>%
  select(-c(X, auc, auc_uncorrected))

# sum_i1_prev <- i1 %>%
#   group_by(prev) %>%
#   summarize(auc_change = mean(auc_diff))
#
# sum_i1_int <- i1 %>%
#   group_by(b_int) %>%
#   summarize(auc_change = mean(auc_diff))

i1_lb <- i1_lb %>%
  filter((auc != 0) & (auc_uncorrected != 0)) %>%
  select(-c(X, auc, auc_uncorrected))

# sum_i1_lb_prev <- i1_lb %>%
#   group_by(prev) %>%
#   summarize(auc_change = mean(auc_diff))
#
# sum_i1_lb_int <- i1_lb %>%
#   group_by(b_int) %>%
#   summarize(auc_change = mean(auc_diff))

i1_df_plot <- i1 %>%
  inner_join(i1_lb, by = c("prev", "b_int", "cmi"), suffix = c(".weighted_auc", ".adjusted_risk")) %>%
  pivot_longer(cols = starts_with("auc"),
               names_to = "method",
               values_to = "auc_change")
i1_df_plot

i1_prev_plot <- ggplot(i1_df_plot, aes(x = method, y = auc_change)) +
  geom_boxplot() +
  facet_grid(.~prev) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "AUC improvement", x = "Adjustment method", y = "AUC change")

i1_prev_plot

i1_int_plot <- ggplot(i1_df_plot, aes(x = method, y = auc_change)) +
  geom_boxplot() +
  facet_grid(.~b_int) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "AUC improvement", x = "Adjustment method", y = "AUC change")

i1_int_plot

# intervention 2
i2 <- read.csv("model_sum2_0216.csv")
i2_lb <- read.csv("model_sum_lb2_0216.csv")

i2 <- i2 %>%
  filter(auc != 0) %>%
  select(-c(X, auc, auc_uncorrected))

# sum_i2_prev <- i2 %>%
#   group_by(prev) %>%
#   summarize(auc_change = mean(auc_diff))
#
# sum_i2_int <- i2 %>%
#   group_by(b_int) %>%
#   summarize(auc_change = mean(auc_diff))

i2_lb <- i2_lb %>%
  filter((auc != 0) & (auc_uncorrected != 0)) %>%
  select(-c(X, auc, auc_uncorrected))

# sum_i1_lb_prev <- i1_lb %>%
#   group_by(prev) %>%
#   summarize(auc_change = mean(auc_diff))
#
# sum_i1_lb_int <- i1_lb %>%
#   group_by(b_int) %>%
#   summarize(auc_change = mean(auc_diff))

i2_df_plot <- i2 %>%
  full_join(i2_lb, by = c("prev", "b_int"), suffix = c(".weighted_auc", ".adjusted_risk")) %>%
  pivot_longer(cols = starts_with("auc"),
               names_to = "method",
               values_to = "auc_change")
i2_df_plot

i2_prev_plot <- ggplot(i2_df_plot, aes(x = method, y = auc_change)) +
  geom_boxplot() +
  facet_grid(.~prev) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "AUC improvement", x = "Adjustment method", y = "AUC change")

i2_prev_plot

i1_int_plot <- ggplot(i1_df_plot, aes(x = method, y = auc_change)) +
  geom_boxplot() +
  facet_grid(.~b_int) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "AUC improvement", x = "Adjustment method", y = "AUC change")

i1_int_plot


