library(foreach)

foreach(a = 1:3, b = rep(10, 3), .combine = '+') %do% rnorm(2)

AUC_list <- list()
MSE_list <- list()
ideal_auc <- c(NA,3)

foreach(p=1:3) %do% {
  AUC.grid <- rep(NA, 10)
  MSE.grid <- rep(NA, 10)
  for (i in 1:10) {
    AUC_i <- 0
    MSE_i <- 0
    for (j in 1:repeats){
      c <- as.data.frame(gen_outcome_f1(risk_t1,
                                        n,
                                        b_intervention = i/10,
                                        cmi_cutoff,
                                        prevalence = p/3))
      ideal_auc[p] <- auc(c$outcome_noint, c$risk_t1)
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
      fit <- glm(y.train ~ risk_t1, data = as.data.frame(x.train), family = binomial(link = "logit"))
      fit_sum <- summary(fit)
      y.pred <- predict(fit, as.data.frame(x.test), type = "response")
      AUC_i = AUC_i + auc(y.test, y.pred)
      MSE_i = MSE_i + mean(fit$residuals^2)
    }
    # AUC
    AUC.grid[i] <- AUC_i/repeats
    MSE.grid[i] <- MSE_i/repeats
  }
  AUC_list[[p]] <- AUC.grid
  MSE_list[[p]] <- MSE.grid
}


for (i in 1:3) {
    for (j in 1:3){

    }
  }
