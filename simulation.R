#============================================================================#
# Code for rare outcome simulation (I + II + III)                            #
# Author: Jiayi Ji                                                           #   
# Date: 07/11/2020                                                           #  
#============================================================================#
library(dplyr)
library(twang)
library(BART)
source("code/bart.multiTrt.R")
source("code/simulation_1_design_rare.R")
source("code/simulation_2_design_rare.R")
source("code/simulation_3_design_rare.R")
source("code/ate_fun.R")
source("code/RAMS.R")
source("code/postSumm.R")
source("code/bart.multiTrt.R")
# Simulation 1: Ratio of units = 1:1:1 n=1500 --------------------------------------------------------------
nsim = 200
set.seed(2)
ate_result_scenario_1 <- matrix(NA, nrow = nsim, ncol = 24)
for (i in 1:nsim) {
  idata <- data_gen(n = 1500, ratio = 1)
  y <- idata$Yobs
  trt_ind<- as.factor(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  psdat<- cbind(trt_ind, all_vars)
  # GBM to get ps
  psmod<-twang::mnps(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
              bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
              verbose = TRUE, estimand = "ATE", stop.method = "es.max",  sampw = NULL,
              treatATT = NULL)
  wt_hat<- get.weights(psmod, stop.method = "es.max",estimand = "ATE")
  ate_gbm_ps1 <- psmod$psList$`1`$ps %>% pull(es.max.ATE)
  ate_gbm_ps2 <- psmod$psList$`2`$ps %>% pull(es.max.ATE)
  # MLR to get ps
  psmod2 <-  nnet::multinom(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = psdat,trace = FALSE)
  pred.ps <- fitted(psmod2)
  ate.ps1 <- pred.ps[,1]
  ate.ps2 <- pred.ps[,2]
  ate.wt.1<- 1/pred.ps[,1]
  ate.wt.2<- 1/pred.ps[,2]
  ate.wt.3<- 1/pred.ps[,3]  
  #1a. to compute ATEs using MLR estimated weights 
  RR12.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[1]]
  RR13.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[2]]
  RR23.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[3]]
  #1b. to compute ATEs using GBM estimated weights
  RR12.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[1]]
  RR13.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[2]]
  RR23.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[3]]
  #2a. to compute ATEs using trimmed MLR weights 
  ate.wt.1.trunc<- trunc_fun(ate.wt.1)
  ate.wt.2.trunc<- trunc_fun(ate.wt.2)
  ate.wt.3.trunc<- trunc_fun(ate.wt.3)
  RR12.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[1]]
  RR13.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[2]]
  RR23.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[3]]
  #2b. to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RR12.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[1]]
  RR13.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[2]]
  RR23.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[3]]
  
  #3a. to compute ATEs spline of GPS using MLR estimated weights
  RR.gps.3a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1, ps2 = ate.ps2)
  
  #3b. to compute ATEs spline of GPS using GBM estimated weights
  RR.gps.3b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm_ps1, ps2 = ate_gbm_ps2)
  #4a. to compute ATEs spline of GPS using trimmed MLR weights
  ate.ps1.trunc <- trunc_fun(ate.ps1)
  ate.ps2.trunc <- trunc_fun(ate.ps2)
  RR.gps.4a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1.trunc, ps2 = ate.ps2.trunc)
  #4b. to compute ATEs spline of GPS using trimmed GBM weights
  
  ate_gbm.ps1.trunc <- trunc_fun(ate_gbm_ps1)
  ate_gbm.ps2.trunc <- trunc_fun(ate_gbm_ps2)
  RR.gps.4b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm.ps1.trunc, ps2 = ate_gbm.ps2.trunc)
  
  ate_result_scenario_1[i,] = c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.1b, RR13.iptw.1b, RR23.iptw.1b, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a, RR12.iptw.2b, RR13.iptw.2b, RR23.iptw.2b, RR.gps.3a, RR.gps.3b, RR.gps.4a, RR.gps.4b)
  
  # ate_res_mlr[i,] <- c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a)
  print(i)
}
# BART
set.seed(2)
for (i in 1:nsim) {
  idata <- data_gen(n = 1500, ratio = 1)
  y <- idata$Yobs
  trt_ind<- as.numeric(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  ate = bart.multiTrt(
    y = y,
    x = all_vars,
    trt = trt_ind,
    estimand = "ATE",
    k = 2,
    ntree = 100,
    ndpost = 1000,
    nskip = 1000
  )
  ate_bart_result_scenario_1 <- rbind(cbind(ate$ATE12, ate$ATE13, ate$ATE23))
  print(i)
}
# Simulation 1：Ratio of units = 1:4:3 n=4000 ------------------------------------------------------------
nsim = 200
ate_result_scenario_2 <- matrix(NA, nrow = nsim, ncol = 24)
set.seed(2)
for (i in 1:nsim) {
  idata <- data_gen(n = 4000, ratio = 2)
  y <- idata$Yobs
  trt_ind<- as.factor(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  psdat<- cbind(trt_ind, all_vars)
  psmod<-twang::mnps(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
                     bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
                     verbose = TRUE, estimand = "ATE", stop.method = "es.max",  sampw = NULL,   
                     treatATT = NULL)
  # GBM to get ps
  wt_hat<- get.weights(psmod, stop.method = "es.max",estimand = "ATE")
  ate_gbm_ps1 <- psmod$psList$`1`$ps %>% pull(es.max.ATE)
  ate_gbm_ps2 <- psmod$psList$`2`$ps %>% pull(es.max.ATE)
  # MLR to get ps
  psmod2 <-  nnet::multinom(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = psdat,trace = FALSE)
  pred.ps <- fitted(psmod2)
  ate.ps1 <- pred.ps[,1]
  ate.ps2 <- pred.ps[,2]
  ate.wt.1<- 1/pred.ps[,1]
  ate.wt.2<- 1/pred.ps[,2]
  ate.wt.3<- 1/pred.ps[,3]  
  #1a. to compute ATEs using MLR estimated weights 
  RR12.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[4]]
  RR13.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[5]]
  RR23.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[6]]
  #1b. to compute ATEs using GBM estimated weights
  RR12.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[4]]
  RR13.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[5]]
  RR23.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[6]]
  #2a. to compute ATEs using trimmed MLR weights 
  ate.wt.1.trunc<- trunc_fun(ate.wt.1)
  ate.wt.2.trunc<- trunc_fun(ate.wt.2)
  ate.wt.3.trunc<- trunc_fun(ate.wt.3)
  RR12.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[4]]
  RR13.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[5]]
  RR23.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[6]]
  #2b. to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RR12.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[4]]
  RR13.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[5]]
  RR23.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[6]]
  
  #3a. to compute ATEs spline of GPS using MLR estimated weights
  RR.gps.3a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1, ps2 = ate.ps2)
  
  #3b. to compute ATEs spline of GPS using GBM estimated weights
  RR.gps.3b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm_ps1, ps2 = ate_gbm_ps2)
  #4a. to compute ATEs spline of GPS using trimmed MLR weights
  ate.ps1.trunc <- trunc_fun(ate.ps1)
  ate.ps2.trunc <- trunc_fun(ate.ps2)
  RR.gps.4a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1.trunc, ps2 = ate.ps2.trunc)
  #4b. to compute ATEs spline of GPS using trimmed GBM weights
  
  ate_gbm.ps1.trunc <- trunc_fun(ate_gbm_ps1)
  ate_gbm.ps2.trunc <- trunc_fun(ate_gbm_ps2)
  RR.gps.4b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm.ps1.trunc, ps2 = ate_gbm.ps2.trunc)
  
  ate_result_scenario_2[i,] = c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.1b, RR13.iptw.1b, RR23.iptw.1b, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a, RR12.iptw.2b, RR13.iptw.2b, RR23.iptw.2b, RR.gps.3a, RR.gps.3b, RR.gps.4a, RR.gps.4b)
  
  # ate_res_mlr[i,] <- c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a)
  print(i)
}

# BART
set.seed(2)
for (i in 1:nsim) {
  idata <- data_gen(n = 4000, ratio = 2)
  y <- idata$Yobs
  trt_ind<- as.numeric(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  ate = bart.multiTrt(
    y = y,
    x = all_vars,
    trt = trt_ind,
    estimand = "ATE",
    k = 2,
    ntree = 100,
    ndpost = 1000,
    nskip = 1000
  )
  ate_bart_result_scenario_1 <- rbind(cbind(ate$ATE12, ate$ATE13, ate$ATE23))
  print(i)
}

# Simulation 1：Ratio of units = 1:10:8 n=9500 ------------------------------------------------------------
nsim = 200
ate_result_scenario_3 <- matrix(NA, nrow = nsim, ncol = 24)
set.seed(2)
for (i in 1:200) {
  idata <- data_gen(n = 9500, ratio = 3)
  y <- idata$Yobs
  trt_ind<- as.factor(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  psdat<- cbind(trt_ind, all_vars)
  # GBM to get ps
  psmod<-twang::mnps(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
                     bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
                     verbose = TRUE, estimand = "ATE", stop.method = "es.max",  sampw = NULL,   
                     treatATT = NULL)
  wt_hat<- get.weights(psmod, stop.method = "es.max",estimand = "ATE")
  ate_gbm_ps1 <- psmod$psList$`1`$ps %>% pull(es.max.ATE)
  ate_gbm_ps2 <- psmod$psList$`2`$ps %>% pull(es.max.ATE)
  # MLR to get ps
  psmod2 <-  nnet::multinom(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = psdat,trace = FALSE)
  pred.ps <- fitted(psmod2)
  ate.ps1 <- pred.ps[,1]
  ate.ps2 <- pred.ps[,2]
  ate.wt.1<- 1/pred.ps[,1]
  ate.wt.2<- 1/pred.ps[,2]
  ate.wt.3<- 1/pred.ps[,3]  
  #1a. to compute ATEs using MLR estimated weights 
  RR12.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[4]]
  RR13.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[5]]
  RR23.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[6]]
  #1b. to compute ATEs using GBM estimated weights
  RR12.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[4]]
  RR13.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[5]]
  RR23.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[6]]
  #2a. to compute ATEs using trimmed MLR weights 
  ate.wt.1.trunc<- trunc_fun(ate.wt.1)
  ate.wt.2.trunc<- trunc_fun(ate.wt.2)
  ate.wt.3.trunc<- trunc_fun(ate.wt.3)
  RR12.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[4]]
  RR13.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[5]]
  RR23.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[6]]
  #2b. to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RR12.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[4]]
  RR13.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[5]]
  RR23.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[6]]
  
  #3a. to compute ATEs spline of GPS using MLR estimated weights
  RR.gps.3a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1, ps2 = ate.ps2)
  
  #3b. to compute ATEs spline of GPS using GBM estimated weights
  RR.gps.3b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm_ps1, ps2 = ate_gbm_ps2)
  #4a. to compute ATEs spline of GPS using trimmed MLR weights
  ate.ps1.trunc <- trunc_fun(ate.ps1)
  ate.ps2.trunc <- trunc_fun(ate.ps2)
  RR.gps.4a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1.trunc, ps2 = ate.ps2.trunc)
  #4b. to compute ATEs spline of GPS using trimmed GBM weights
  
  ate_gbm.ps1.trunc <- trunc_fun(ate_gbm_ps1)
  ate_gbm.ps2.trunc <- trunc_fun(ate_gbm_ps2)
  RR.gps.4b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm.ps1.trunc, ps2 = ate_gbm.ps2.trunc)
  
  ate_result_scenario_3[i,] = c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.1b, RR13.iptw.1b, RR23.iptw.1b, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a, RR12.iptw.2b, RR13.iptw.2b, RR23.iptw.2b, RR.gps.3a, RR.gps.3b, RR.gps.4a, RR.gps.4b)
  
  # ate_res_mlr[i,] <- c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a)
  print(i)
}

# BART
set.seed(2)
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, ratio = 3)
  y <- idata$Yobs
  trt_ind<- as.numeric(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  ate = bart.multiTrt(
    y = y,
    x = all_vars,
    trt = trt_ind,
    estimand = "ATE",
    k = 2,
    ntree = 100,
    ndpost = 1000,
    nskip = 1000
  )
  ate_bart_result_scenario_1 <- rbind(cbind(ate$ATE12, ate$ATE13, ate$ATE23))
  print(i)
}

# Simulation 2：Outcome prevalence = 1-5%  n=9500------------------------------------------------------------
nsim = 200
set.seed(2)
ate_result_scenario_1 <- matrix(NA, nrow = nsim, ncol = 24)
# strFormula  = "trt_ind~X1    +      X2   +      X3  +       X4     +     X5+ X6+ X7+ X8+ X9 +X10"

for (i in 1:nsim) {
  idata <- data_gen(n = 9500, outc.prev = "1%-5%")
  y <- idata$Yobs
  trt_ind<- as.factor(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  psdat<- cbind(trt_ind, all_vars)
  # GBM to get ps
  psmod<-twang::mnps(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
                     bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
                     verbose = TRUE, estimand = "ATE", stop.method = "es.max",  sampw = NULL,
                     treatATT = NULL)
  wt_hat<- get.weights(psmod, stop.method = "es.max",estimand = "ATE")
  ate_gbm_ps1 <- psmod$psList$`1`$ps %>% pull(es.max.ATE)
  ate_gbm_ps2 <- psmod$psList$`2`$ps %>% pull(es.max.ATE)
  # MLR to get ps
  psmod2 <-  nnet::multinom(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = psdat,trace = FALSE)
  pred.ps <- fitted(psmod2)
  ate.ps1 <- pred.ps[,1]
  ate.ps2 <- pred.ps[,2]
  ate.wt.1<- 1/pred.ps[,1]
  ate.wt.2<- 1/pred.ps[,2]
  ate.wt.3<- 1/pred.ps[,3]  
  #1a. to compute ATEs using MLR estimated weights 
  RR12.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[4]]
  RR13.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[5]]
  RR23.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[6]]
  #1b. to compute ATEs using GBM estimated weights
  RR12.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[4]]
  RR13.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[5]]
  RR23.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[6]]
  #2a. to compute ATEs using trimmed MLR weights 
  ate.wt.1.trunc<- trunc_fun(ate.wt.1)
  ate.wt.2.trunc<- trunc_fun(ate.wt.2)
  ate.wt.3.trunc<- trunc_fun(ate.wt.3)
  RR12.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[4]]
  RR13.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[5]]
  RR23.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[6]]
  #2b. to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RR12.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[4]]
  RR13.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[5]]
  RR23.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[6]]
  
  #3a. to compute ATEs spline of GPS using MLR estimated weights
  RR.gps.3a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1, ps2 = ate.ps2)
  
  #3b. to compute ATEs spline of GPS using GBM estimated weights
  RR.gps.3b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm_ps1, ps2 = ate_gbm_ps2)
  #4a. to compute ATEs spline of GPS using trimmed MLR weights
  ate.ps1.trunc <- trunc_fun(ate.ps1)
  ate.ps2.trunc <- trunc_fun(ate.ps2)
  RR.gps.4a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1.trunc, ps2 = ate.ps2.trunc)
  #4b. to compute ATEs spline of GPS using trimmed GBM weights
  
  ate_gbm.ps1.trunc <- trunc_fun(ate_gbm_ps1)
  ate_gbm.ps2.trunc <- trunc_fun(ate_gbm_ps2)
  RR.gps.4b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm.ps1.trunc, ps2 = ate_gbm.ps2.trunc)
  
  ate_result_scenario_1[i,] = c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.1b, RR13.iptw.1b, RR23.iptw.1b, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a, RR12.iptw.2b, RR13.iptw.2b, RR23.iptw.2b, RR.gps.3a, RR.gps.3b, RR.gps.4a, RR.gps.4b)
  
  # ate_res_mlr[i,] <- c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a)
  print(i)
}
# BART
nsim = 200
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, outc.prev = "1%-5%")
  y <- idata$Yobs
  trt_ind<- as.numeric(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  ate = bart.multiTrt(
    y = y,
    x = all_vars,
    trt = trt_ind,
    estimand = "ATE",
    k = 2,
    ntree = 100,
    ndpost = 1000,
    nskip = 1000
  )
  ate_bart_result_scenario_2 <- rbind(cbind(ate$ATE12, ate$ATE13, ate$ATE23))
  save(ate_bart_result_scenario_2, file= paste0("Rdata/bart_sim_rare_p2_scenario_2_no_discard_",i,".RData"))
  print(i)
}
# Simulation 2：Outcome prevalence = 5-10%  n=9500------------------------------------------------------------
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, outc.prev = "5%-10%")
  y <- idata$Yobs
  trt_ind<- as.factor(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  psdat<- cbind(trt_ind, all_vars)
  # GBM to get ps
  psmod<-twang::mnps(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
                     bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
                     verbose = TRUE, estimand = "ATE", stop.method = "es.max",  sampw = NULL,
                     treatATT = NULL)
  wt_hat<- get.weights(psmod, stop.method = "es.max",estimand = "ATE")
  ate_gbm_ps1 <- psmod$psList$`1`$ps %>% pull(es.max.ATE)
  ate_gbm_ps2 <- psmod$psList$`2`$ps %>% pull(es.max.ATE)
  # MLR to get ps
  psmod2 <-  nnet::multinom(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = psdat,trace = FALSE)
  pred.ps <- fitted(psmod2)
  ate.ps1 <- pred.ps[,1]
  ate.ps2 <- pred.ps[,2]
  ate.wt.1<- 1/pred.ps[,1]
  ate.wt.2<- 1/pred.ps[,2]
  ate.wt.3<- 1/pred.ps[,3]  
  #1a. to compute ATEs using MLR estimated weights 
  RR12.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[4]]
  RR13.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[5]]
  RR23.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[6]]
  #1b. to compute ATEs using GBM estimated weights
  RR12.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[4]]
  RR13.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[5]]
  RR23.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[6]]
  #2a. to compute ATEs using trimmed MLR weights 
  ate.wt.1.trunc<- trunc_fun(ate.wt.1)
  ate.wt.2.trunc<- trunc_fun(ate.wt.2)
  ate.wt.3.trunc<- trunc_fun(ate.wt.3)
  RR12.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[4]]
  RR13.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[5]]
  RR23.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[6]]
  #2b. to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RR12.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[4]]
  RR13.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[5]]
  RR23.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[6]]
  
  #3a. to compute ATEs spline of GPS using MLR estimated weights
  RR.gps.3a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1, ps2 = ate.ps2)
  
  #3b. to compute ATEs spline of GPS using GBM estimated weights
  RR.gps.3b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm_ps1, ps2 = ate_gbm_ps2)
  #4a. to compute ATEs spline of GPS using trimmed MLR weights
  ate.ps1.trunc <- trunc_fun(ate.ps1)
  ate.ps2.trunc <- trunc_fun(ate.ps2)
  RR.gps.4a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1.trunc, ps2 = ate.ps2.trunc)
  #4b. to compute ATEs spline of GPS using trimmed GBM weights
  
  ate_gbm.ps1.trunc <- trunc_fun(ate_gbm_ps1)
  ate_gbm.ps2.trunc <- trunc_fun(ate_gbm_ps2)
  RR.gps.4b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm.ps1.trunc, ps2 = ate_gbm.ps2.trunc)
  
  ate_result_scenario_3[i,] = c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.1b, RR13.iptw.1b, RR23.iptw.1b, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a, RR12.iptw.2b, RR13.iptw.2b, RR23.iptw.2b, RR.gps.3a, RR.gps.3b, RR.gps.4a, RR.gps.4b)
  
  # ate_res_mlr[i,] <- c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a)
  print(i)
}
# BART
nsim = 200
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, outc.prev = "5%-10%")
  y <- idata$Yobs
  trt_ind<- as.numeric(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  ate = bart.multiTrt(
    y = y,
    x = all_vars,
    trt = trt_ind,
    estimand = "ATE",
    k = 2,
    ntree = 100,
    ndpost = 1000,
    nskip = 1000
  )
  ate_bart_result_scenario_2 <- rbind(cbind(ate$ATE12, ate$ATE13, ate$ATE23))
  save(ate_bart_result_scenario_2, file= paste0("Rdata/bart_sim_rare_p2_scenario_2_no_discard_",i,".RData"))
  print(i)
}

# Simulation 3：Strong overlap n=9500------------------------------------------------------------
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, overlap = "extremely strong")
  y <- idata$Yobs
  trt_ind<- as.factor(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  psdat<- cbind(trt_ind, all_vars)
  # GBM to get ps
  psmod<-twang::mnps(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
                     bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
                     verbose = TRUE, estimand = "ATE", stop.method = "es.max",  sampw = NULL,
                     treatATT = NULL)
  wt_hat<- get.weights(psmod, stop.method = "es.max",estimand = "ATE")
  ate_gbm_ps1 <- psmod$psList$`1`$ps %>% pull(es.max.ATE)
  ate_gbm_ps2 <- psmod$psList$`2`$ps %>% pull(es.max.ATE)
  # MLR to get ps
  psmod2 <-  nnet::multinom(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = psdat,trace = FALSE)
  pred.ps <- fitted(psmod2)
  ate.ps1 <- pred.ps[,1]
  ate.ps2 <- pred.ps[,2]
  ate.wt.1<- 1/pred.ps[,1]
  ate.wt.2<- 1/pred.ps[,2]
  ate.wt.3<- 1/pred.ps[,3]  
  #1a. to compute ATEs using MLR estimated weights 
  RR12.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[4]]
  RR13.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[5]]
  RR23.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[6]]
  #1b. to compute ATEs using GBM estimated weights
  RR12.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[4]]
  RR13.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[5]]
  RR23.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[6]]
  #2a. to compute ATEs using trimmed MLR weights 
  ate.wt.1.trunc<- trunc_fun(ate.wt.1)
  ate.wt.2.trunc<- trunc_fun(ate.wt.2)
  ate.wt.3.trunc<- trunc_fun(ate.wt.3)
  RR12.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[4]]
  RR13.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[5]]
  RR23.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[6]]
  #2b. to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RR12.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[4]]
  RR13.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[5]]
  RR23.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[6]]
  
  #3a. to compute ATEs spline of GPS using MLR estimated weights
  RR.gps.3a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1, ps2 = ate.ps2)
  
  #3b. to compute ATEs spline of GPS using GBM estimated weights
  RR.gps.3b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm_ps1, ps2 = ate_gbm_ps2)
  #4a. to compute ATEs spline of GPS using trimmed MLR weights
  ate.ps1.trunc <- trunc_fun(ate.ps1)
  ate.ps2.trunc <- trunc_fun(ate.ps2)
  RR.gps.4a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1.trunc, ps2 = ate.ps2.trunc)
  #4b. to compute ATEs spline of GPS using trimmed GBM weights
  
  ate_gbm.ps1.trunc <- trunc_fun(ate_gbm_ps1)
  ate_gbm.ps2.trunc <- trunc_fun(ate_gbm_ps2)
  RR.gps.4b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm.ps1.trunc, ps2 = ate_gbm.ps2.trunc)
  
  ate_result_scenario_2[i,] = c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.1b, RR13.iptw.1b, RR23.iptw.1b, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a, RR12.iptw.2b, RR13.iptw.2b, RR23.iptw.2b, RR.gps.3a, RR.gps.3b, RR.gps.4a, RR.gps.4b)
  
  # ate_res_mlr[i,] <- c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a)
  print(i)
}
nsim = 200
set.seed(2)
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, overlap = "extremely strong")
  y <- idata$Yobs
  trt_ind<- as.numeric(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  ate = bart.multiTrt(
    y = y,
    x = all_vars,
    trt = trt_ind,
    estimand = "ATE",
    k = 2,
    ntree = 100,
    ndpost = 1000,
    nskip = 1000
  )
  ate_bart_result_scenario_1 <- rbind(cbind(ate$ATE12, ate$ATE13, ate$ATE23))
  print(i)
}

# Simulation 3：Moderate overlap n=9500------------------------------------------------------------
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, overlap = "strong")
  y <- idata$Yobs
  trt_ind<- as.factor(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  psdat<- cbind(trt_ind, all_vars)
  # GBM to get ps
  psmod<-twang::mnps(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
                     bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
                     verbose = TRUE, estimand = "ATE", stop.method = "es.max",  sampw = NULL,
                     treatATT = NULL)
  wt_hat<- get.weights(psmod, stop.method = "es.max",estimand = "ATE")
  ate_gbm_ps1 <- psmod$psList$`1`$ps %>% pull(es.max.ATE)
  ate_gbm_ps2 <- psmod$psList$`2`$ps %>% pull(es.max.ATE)
  # MLR to get ps
  psmod2 <-  nnet::multinom(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = psdat,trace = FALSE)
  pred.ps <- fitted(psmod2)
  ate.ps1 <- pred.ps[,1]
  ate.ps2 <- pred.ps[,2]
  ate.wt.1<- 1/pred.ps[,1]
  ate.wt.2<- 1/pred.ps[,2]
  ate.wt.3<- 1/pred.ps[,3]  
  #1a. to compute ATEs using MLR estimated weights 
  RR12.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[4]]
  RR13.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[5]]
  RR23.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[6]]
  #1b. to compute ATEs using GBM estimated weights
  RR12.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[4]]
  RR13.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[5]]
  RR23.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[6]]
  #2a. to compute ATEs using trimmed MLR weights 
  ate.wt.1.trunc<- trunc_fun(ate.wt.1)
  ate.wt.2.trunc<- trunc_fun(ate.wt.2)
  ate.wt.3.trunc<- trunc_fun(ate.wt.3)
  RR12.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[4]]
  RR13.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[5]]
  RR23.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[6]]
  #2b. to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RR12.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[4]]
  RR13.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[5]]
  RR23.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[6]]
  
  #3a. to compute ATEs spline of GPS using MLR estimated weights
  RR.gps.3a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1, ps2 = ate.ps2)
  
  #3b. to compute ATEs spline of GPS using GBM estimated weights
  RR.gps.3b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm_ps1, ps2 = ate_gbm_ps2)
  #4a. to compute ATEs spline of GPS using trimmed MLR weights
  ate.ps1.trunc <- trunc_fun(ate.ps1)
  ate.ps2.trunc <- trunc_fun(ate.ps2)
  RR.gps.4a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1.trunc, ps2 = ate.ps2.trunc)
  #4b. to compute ATEs spline of GPS using trimmed GBM weights
  
  ate_gbm.ps1.trunc <- trunc_fun(ate_gbm_ps1)
  ate_gbm.ps2.trunc <- trunc_fun(ate_gbm_ps2)
  RR.gps.4b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm.ps1.trunc, ps2 = ate_gbm.ps2.trunc)
  
  ate_result_scenario_2[i,] = c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.1b, RR13.iptw.1b, RR23.iptw.1b, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a, RR12.iptw.2b, RR13.iptw.2b, RR23.iptw.2b, RR.gps.3a, RR.gps.3b, RR.gps.4a, RR.gps.4b)
  
  # ate_res_mlr[i,] <- c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a)
  print(i)
}

nsim = 200
set.seed(2)
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, overlap = "strong")
  y <- idata$Yobs
  trt_ind<- as.numeric(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  ate = bart.multiTrt(
    y = y,
    x = all_vars,
    trt = trt_ind,
    estimand = "ATE",
    k = 2,
    ntree = 100,
    ndpost = 1000,
    nskip = 1000
  )
  ate_bart_result_scenario_1 <- rbind(cbind(ate$ATE12, ate$ATE13, ate$ATE23))
  print(i)
}
# Simulation 3：Weak overlap n=9500------------------------------------------------------------
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, overlap = "moderate")
  y <- idata$Yobs
  trt_ind<- as.factor(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  psdat<- cbind(trt_ind, all_vars)
  # GBM to get ps
  psmod<-twang::mnps(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
                     bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
                     verbose = TRUE, estimand = "ATE", stop.method = "es.max",  sampw = NULL,
                     treatATT = NULL)
  wt_hat<- get.weights(psmod, stop.method = "es.max",estimand = "ATE")
  ate_gbm_ps1 <- psmod$psList$`1`$ps %>% pull(es.max.ATE)
  ate_gbm_ps2 <- psmod$psList$`2`$ps %>% pull(es.max.ATE)
  # MLR to get ps
  psmod2 <-  nnet::multinom(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = psdat,trace = FALSE)
  pred.ps <- fitted(psmod2)
  ate.ps1 <- pred.ps[,1]
  ate.ps2 <- pred.ps[,2]
  ate.wt.1<- 1/pred.ps[,1]
  ate.wt.2<- 1/pred.ps[,2]
  ate.wt.3<- 1/pred.ps[,3]  
  #1a. to compute ATEs using MLR estimated weights 
  RR12.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[4]]
  RR13.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[5]]
  RR23.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[6]]
  #1b. to compute ATEs using GBM estimated weights
  RR12.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[4]]
  RR13.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[5]]
  RR23.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[6]]
  #2a. to compute ATEs using trimmed MLR weights 
  ate.wt.1.trunc<- trunc_fun(ate.wt.1)
  ate.wt.2.trunc<- trunc_fun(ate.wt.2)
  ate.wt.3.trunc<- trunc_fun(ate.wt.3)
  RR12.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[4]]
  RR13.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[5]]
  RR23.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[6]]
  #2b. to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RR12.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[4]]
  RR13.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[5]]
  RR23.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[6]]
  
  #3a. to compute ATEs spline of GPS using MLR estimated weights
  RR.gps.3a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1, ps2 = ate.ps2)
  
  #3b. to compute ATEs spline of GPS using GBM estimated weights
  RR.gps.3b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm_ps1, ps2 = ate_gbm_ps2)
  #4a. to compute ATEs spline of GPS using trimmed MLR weights
  ate.ps1.trunc <- trunc_fun(ate.ps1)
  ate.ps2.trunc <- trunc_fun(ate.ps2)
  RR.gps.4a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1.trunc, ps2 = ate.ps2.trunc)
  #4b. to compute ATEs spline of GPS using trimmed GBM weights
  
  ate_gbm.ps1.trunc <- trunc_fun(ate_gbm_ps1)
  ate_gbm.ps2.trunc <- trunc_fun(ate_gbm_ps2)
  RR.gps.4b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm.ps1.trunc, ps2 = ate_gbm.ps2.trunc)
  
  ate_result_scenario_2[i,] = c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.1b, RR13.iptw.1b, RR23.iptw.1b, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a, RR12.iptw.2b, RR13.iptw.2b, RR23.iptw.2b, RR.gps.3a, RR.gps.3b, RR.gps.4a, RR.gps.4b)
  
  # ate_res_mlr[i,] <- c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a)
  print(i)
}

nsim = 200
set.seed(2)
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, overlap = "moderate")
  y <- idata$Yobs
  trt_ind<- as.numeric(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  ate = bart.multiTrt(
    y = y,
    x = all_vars,
    trt = trt_ind,
    estimand = "ATE",
    k = 2,
    ntree = 100,
    ndpost = 1000,
    nskip = 1000
  )
  ate_bart_result_scenario_1 <- rbind(cbind(ate$ATE12, ate$ATE13, ate$ATE23))
  print(i)
}
# Simulation 3：Near zero overlap n=9500------------------------------------------------------------
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, overlap = "weak")
  y <- idata$Yobs
  trt_ind<- as.factor(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  psdat<- cbind(trt_ind, all_vars)
  # GBM to get ps
  psmod<-twang::mnps(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
                     bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
                     verbose = TRUE, estimand = "ATE", stop.method = "es.max",  sampw = NULL,
                     treatATT = NULL)
  wt_hat<- get.weights(psmod, stop.method = "es.max",estimand = "ATE")
  ate_gbm_ps1 <- psmod$psList$`1`$ps %>% pull(es.max.ATE)
  ate_gbm_ps2 <- psmod$psList$`2`$ps %>% pull(es.max.ATE)
  # MLR to get ps
  psmod2 <-  nnet::multinom(trt_ind ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = psdat,trace = FALSE)
  pred.ps <- fitted(psmod2)
  ate.ps1 <- pred.ps[,1]
  ate.ps2 <- pred.ps[,2]
  ate.wt.1<- 1/pred.ps[,1]
  ate.wt.2<- 1/pred.ps[,2]
  ate.wt.3<- 1/pred.ps[,3]  
  #1a. to compute ATEs using MLR estimated weights 
  RR12.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[4]]
  RR13.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[5]]
  RR23.iptw.1a = ate_fun(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3)[[6]]
  #1b. to compute ATEs using GBM estimated weights
  RR12.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[4]]
  RR13.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[5]]
  RR23.iptw.1b = ate_fun(wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[6]]
  #2a. to compute ATEs using trimmed MLR weights 
  ate.wt.1.trunc<- trunc_fun(ate.wt.1)
  ate.wt.2.trunc<- trunc_fun(ate.wt.2)
  ate.wt.3.trunc<- trunc_fun(ate.wt.3)
  RR12.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[4]]
  RR13.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[5]]
  RR23.iptw.2a = ate_fun(wt1 = ate.wt.1.trunc, wt2 = ate.wt.2.trunc, wt3 = ate.wt.3.trunc)[[6]]
  #2b. to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RR12.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[4]]
  RR13.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[5]]
  RR23.iptw.2b = ate_fun(wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[6]]
  
  #3a. to compute ATEs spline of GPS using MLR estimated weights
  RR.gps.3a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1, ps2 = ate.ps2)
  
  #3b. to compute ATEs spline of GPS using GBM estimated weights
  RR.gps.3b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm_ps1, ps2 = ate_gbm_ps2)
  #4a. to compute ATEs spline of GPS using trimmed MLR weights
  ate.ps1.trunc <- trunc_fun(ate.ps1)
  ate.ps2.trunc <- trunc_fun(ate.ps2)
  RR.gps.4a <- RegSpline(y = y, trt = trt_ind, ps1 = ate.ps1.trunc, ps2 = ate.ps2.trunc)
  #4b. to compute ATEs spline of GPS using trimmed GBM weights
  
  ate_gbm.ps1.trunc <- trunc_fun(ate_gbm_ps1)
  ate_gbm.ps2.trunc <- trunc_fun(ate_gbm_ps2)
  RR.gps.4b <- RegSpline(y = y, trt = trt_ind, ps1 = ate_gbm.ps1.trunc, ps2 = ate_gbm.ps2.trunc)
  
  ate_result_scenario_2[i,] = c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.1b, RR13.iptw.1b, RR23.iptw.1b, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a, RR12.iptw.2b, RR13.iptw.2b, RR23.iptw.2b, RR.gps.3a, RR.gps.3b, RR.gps.4a, RR.gps.4b)
  
  # ate_res_mlr[i,] <- c(RR12.iptw.1a, RR13.iptw.1a, RR23.iptw.1a, RR12.iptw.2a, RR13.iptw.2a, RR23.iptw.2a)
  print(i)
}

nsim = 200
set.seed(2)
for (i in 1:nsim) {
  idata <- data_gen(n = 9500, overlap = "weak")
  y <- idata$Yobs
  trt_ind<- as.numeric(idata$trtdat$trt_ind)
  all_vars <- idata$trtdat[,-1] #exclude treatment indicator
  ate = bart.multiTrt(
    y = y,
    x = all_vars,
    trt = trt_ind,
    estimand = "ATE",
    k = 2,
    ntree = 100,
    ndpost = 1000,
    nskip = 1000
  )
  ate_bart_result_scenario_1 <- rbind(cbind(ate$ATE12, ate$ATE13, ate$ATE23))
  print(i)
}