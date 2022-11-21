library(tidyverse)
library(dplyr)
library(magrittr)
library(ggplot2)
library(prodlim)
library(pROC)
library(predtools)
library(boot)
library("lattice")

# set up -----------------------------------------------------------------------
set.seed(1)
n <- 1000
repeats <- 20

# generate risk_t1
risk_t1 <- runif(n, 0, 1)

# intervention
cmi_cutoff <- 0.5

# intervention effectiveness/probability of the intervention works (fixed, varies in different trials)
# b_intervention <- c(seq(0, 1, 0.1)) # set to 1  for null case
b_intervention <- 0.1

gen_outcome1 <- function(risk_t1, n = n){
  outcome_noint <- rbinom(n, 1, 0.9 * risk_t1) # the constant vary to control true model AUC
  df <- as.data.frame(cbind(risk_t1, outcome_noint))
  return (df)
}

a <- gen_outcome1(risk_t1, n)
auc(a$outcome_noint, a$risk_t1)

# generate cmi (binary)--------------------------------
gen_cmi <- function(risk_t1, n = n, cmi_cutoff){
  # 1. cutoff approach
  df <- gen_outcome1(risk_t1, n)
  df <- df %>%
    mutate(cmi = ifelse(risk_t1 > cmi_cutoff, 1, 0))

  # 2. percentage approach
  # df <- gen_outcome1(risk_t1, n)
  # cmi <- sample(c(0, 1),size = nrow(df), prob = c(1 - cmi_perc, cmi_perc),
  #               replace = TRUE)
  return (df)
}

# b <- gen_cmi(risk_t1, n, cmi_cutoff)


# generate risk_t2 depending on cmi & risk_t1-----------------------------------
gen_risk_t2 <- function(risk_t1, n = n, b_intervention, cmi_cutoff){
  # if cmi == 1 & outcome_noint == 1, then risk_t2 = risk_t1 * b_intervention
  df <- gen_cmi(risk_t1, n, cmi_cutoff)
  df <- df %>%
    mutate(risk_t2 = ifelse((outcome_noint == 1 & cmi == 1), risk_t1 * b_intervention,
                            risk_t1))
  return (df)
}

# c <- gen_risk_t2(risk_t1, n, b_intervention, cmi_cutoff)

# generate outcome from risk_t2-------------------------------------------------
gen_outcome <- function(risk_t1, n = n, b_intervention,cmi_cutoff){
  df <- gen_risk_t2(risk_t1, n, b_intervention, cmi_cutoff)
  condition <- df %>%
    filter(outcome_noint == 1 & cmi == 1)

  # 1. generate all observed outcome probabilistically
  # oo: if outcome_no_intervention = 0 & intervention = 0 --> old outcome
  # ol: if outcome_no_intervention = 0 & intervention = 1 --> old outcome
  # lo: if outcome_no_intervention = 1 & intervention = 0 --> old outcome
  # ll: if outcome_no_intervention = 1 & intervention = 1 --> rbinom(n, 1, risks_t2)

  num_obs <- dim(condition)[1]

  df_final <- df %>%
    mutate(outcome_int = ifelse((outcome_noint  == 1 & cmi == 1),
                                 rbinom(num_obs, 1, df$risk_t2),
                                 outcome_noint))


  # 2. choose b_intervention% of the ppl (outcome_no_intervention = 1 & intervention = 1) to have
  # outcome = 0

  # choose (1-b_intervention)% as fraction who have effective treatment
  # b_intervention = 0.1 (very effective treatment), then 90% change their 1 to 0
  # effective <- sample_frac(condition, size = 1-b_intervention, replace = FALSE)
  # effective_df <- subset(df,(risk_t1 %in% effective$risk_t1 &
  #              outcome_noint %in% effective$outcome_noint &
  #              cmi %in% effective$cmi &
  #              risk_t2 %in% effective$risk_t2))
  #
  # non_effective_df <- subset(df,!(risk_t1 %in% effective$risk_t1 &
  #                                  outcome_noint %in% effective$outcome_noint &
  #                                  cmi %in% effective$cmi &
  #                                  risk_t2 %in% effective$risk_t2))
  # effective_df$outcome_int <- 0
  # non_effective_df$outcome_int <- non_effective_df$outcome_noint
  # df_final <- rbind(effective_df, non_effective_df)

  return (df_final)
}

# d <- gen_outcome(risk_t1, n, b_intervention, cmi_cutoff)


# model training-----------------------------------------------------------------
c <- as.data.frame(gen_outcome(risk_t1, n = n, b_intervention, cmi_cutoff))
x <- c %>%
  select(-outcome_int)
y <- c %>%
  select(outcome_int)

# train/test split
train <- sample(1:nrow(x), nrow(x)/2)
x.train <- x[train,]
y.train <- y[train,]

test <- (-train)
x.test <- x[test,]
y.test <- y[test,]

# fit logistic regression
fit <- glm(y.train ~ risk_t1, data = as.data.frame(x.train), family = "binomial")

# Evaluate-----------------------------------------------------------------------
y.pred <- predict(fit, as.data.frame(x.test), type = "response")

# AUC
AUC <- auc(y.test, y.pred)
AUC

# calibration
test_data <- as.data.frame(cbind(y.test, y.pred))
calibration_fit <- lm(y.test ~ y.pred, data = test_data)
slope <- calibration_fit$coefficients[2]

ggplot(data = test_data,aes(x = y.pred, y = y.test)) +
  xlim(0,1) +
  ylim(0,1) +
  geom_point()+
  geom_smooth(method = lm, se = FALSE)+
  geom_abline(intercept = 0, slope = 1, color = "red")


# AUC grid----------------------------------------------------------------------
# vary b_intervention
AUC.grid <- rep(NA, 9)
for (i in 1:9) {
  AUC_i <- 0
  for (j in 1:repeats){
  c <- as.data.frame(gen_outcome(risk_t1,
                                 n,
                                 b_intervention = i/10,
                                 cmi_cutoff))
  print(mean(c$outcome_int))
  x <- c %>%
    select(risk_t1)
  y <- c %>%
    select(outcome_int)

  # train/test split
  train <- sample(1:nrow(x), nrow(x)/2)
  x.train <- x[train,]
  y.train <- y[train,]

  test <- (-train)
  x.test <- x[test,]
  y.test <- y[test,]

  fit <- glm(y.train ~ ., data = as.data.frame(x.train), family = "binomial")
  y.pred <- predict(fit, as.data.frame(x.test), type = "response")
  AUC_i = AUC_i + auc(y.test, y.pred)
  }
  # AUC
  AUC.grid[i] <- AUC_i/repeats
}

plot(AUC.grid, xlab = "b_intervention", xaxt = "n", ylab = "AUC", ylim = c(0.5, 1),
     type = "b", sub = paste0("cmi_cutoff =", cmi_cutoff))
axis(side=1, at=1:10, labels = seq(0.1, 1, 0.1))

# vary cmi_cutoff
cmi.grid <- rep(NA, 9)
for (i in 1:9) {
  cmi_i <- 0
  for (j in 1:repeats){
  c <- as.data.frame(gen_outcome(risk_t1,
                                 n,
                                 b_intervention,
                                 cmi_cutoff = i/10))
  x <- c %>%
    select(risk_t1)
  y <- c %>%
    select(outcome_int)

  # train/test split
  train <- sample(1:nrow(x), nrow(x)/2)
  x.train <- x[train,]
  y.train <- y[train,]

  test <- (-train)
  x.test <- x[test,]
  y.test <- y[test,]

  fit <- glm(y.train ~ ., data = as.data.frame(x.train), family = "binomial")
  y.pred <- predict(fit, as.data.frame(x.test), type = "response")
  cmi_i = cmi_i + auc(y.test, y.pred)
  }
  cmi.grid[i] <- cmi_i/repeats
}

plot(cmi.grid, xlab = "CMI_cutoff", xaxt = "n", ylab = "AUC", ylim = c(0.5, 1),
     type = "b", sub = paste0("b_intervention =", b_intervention))
axis(side=1, at=1:10, labels = seq(0.1, 1, 0.1))


# heatmap -----------------------------------------------------------------
cmi <- seq(0.1, 0.9, length.out=9)
b_int <- seq(0.1, 0.9, length.out=9)

data_auc <- expand.grid(X=cmi, Y=b_int)
data_mse <- expand.grid(X=cmi, Y=b_int)
data_calibration <- expand.grid(X=cmi, Y=b_int)

auc_table <- matrix(NA, nrow = length(cmi), ncol = length(cmi))
mse_table <- matrix(NA, nrow = length(cmi), ncol = length(cmi))
calibration_table <- matrix(NA, nrow = length(cmi), ncol = length(cmi))

for (i in 1:9){
  for (j in 1:9){
    cur_auc <- 0
    cur_mse <- 0
    cur_calibration <- 0
    for (k in 1:repeats){
      c <- as.data.frame(gen_outcome(risk_t1,
                                     n,
                                     b_intervention = i/10,
                                     cmi_cutoff = j/10))
      x <- c %>%
        select(-outcome_int)
      y <- c %>%
        select(outcome_int)

      # train/test split
      train <- sample(1:nrow(x), nrow(x)/2)
      x.train <- x[train,]
      y.train <- y[train,]

      test <- (-train)
      x.test <- x[test,]
      y.test <- y[test,]

      fit <- glm(y.train ~ risk_t1, data = x.train, family = "binomial")
      y.pred <- predict(fit, as.data.frame(x.test), type = "response")

      # update auc and mse
      cur_auc = cur_auc + auc(y.test, y.pred)
      cur_mse = cur_mse + sqrt(mean((y.test-y.pred)^2))

      test_data <- as.data.frame(cbind(y.test, y.pred))
      calibration_fit <- lm(y.test ~ y.pred, data = test_data)

      cur_calibration = cur_calibration + ifelse(calibration_fit$coefficients[2] < -600,
                                                 0, calibration_fit$coefficients[2])
    }

    auc_table[i,j] <- cur_auc/repeats
    mse_table[i,j] <- cur_mse/repeats
    calibration_table[i,j] <- cur_calibration/repeats
  }
}

auc_plot <- levelplot(auc_table ~ X*Y, data=data_auc,
                      main="AUC",
                      xlab="cmi", ylab="b_int",
                      col.regions = heat.colors(100))
auc_plot

mse_plot <- levelplot(mse_table ~ X*Y, data=data_mse,
                      main="Mean Squared Error",
                      xlab="cmi", ylab="b_int",
                      col.regions = heat.colors(100))
mse_plot


calibration_plot <- levelplot(calibration_table ~ X*Y, data=data_calibration,
                      main="Calibration curve slope",
                      xlab="cmi", ylab="b_int",
                      col.regions = cm.colors(100))
calibration_plot


