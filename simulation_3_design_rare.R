#============================================================================#
# Simulation design - part 3                                                 #
# sample size 9500 with ratio of units 1:10:8                                #
# Covariate overlap "weak", "moderate I", "moderate II" and "strong"         #
# number of confounders 10                                                   #
# Author: Liangyuan Hu & Chenyang Gu                                         #
# Feb 11, 2020                                                          #
#============================================================================#

expit = function(x) {exp(x)/(1+exp(x))}

data_gen = function(n = 9500, overlap = "weak", all_confounder = T) {

    p = 10

    #generate treatment label W
    W = sample(c(1,2,3), size=n, replace=T, prob=c(.05,.53,.42))

    #==================================#
    #===========weak overlap===========#
    #==================================#
    if (overlap == "weak") {
        Xcon = matrix(NA, nrow=n, ncol=5)
        Xcat = matrix(NA, nrow=n, ncol=5)

        Xcon<-matrix(NA, nrow=n, ncol=5); Xcat<-matrix(NA, nrow=n, ncol=5);
        for (i in 1:n){
          Xcon[i,] = matrix(rnorm(5,-1,1)*(W[i]==1) + rnorm(5, 1, 1)*(W[i]==2) + rnorm(5, 3, 1)*(W[i]==3),ncol=5);
          Xcat[i,] = matrix(sample (c(1,2,3), size = 5,replace =T, prob = c(.3,.3,.4))*(W[i]==1)
                            + sample (c(1,2,3), size = 5,replace =T, prob = c(.6,.2,.2))*(W[i]==2)
                            + sample (c(1,2,3), size = 5,replace =T, prob = c(.8,.1,.1))*(W[i]==3),ncol=5)
        }
        x1 = Xcon[,1]; x2 = Xcon[,2]; x3 = Xcon[,3]; x4 = Xcon[,4]; x5 = Xcon[,5];
        x6 = Xcat[,1]; x7 = Xcat[,2]; x8 = Xcat[,3]; x9 = Xcat[,4]; x10 = Xcat[,5]
        trtdat<- data.frame(W,Xcon, Xcat)


        #parallel response surface model
        #set.seed(3242019)
        tau1 = -6.0; tau2 = -4.0; tau3 = -5.0
        Yp1 = expit(tau1 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)
        Yp2 = expit(tau2 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)
        Yp3 = expit(tau3 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)

        #potential outcomes
        Y1 = Y2 = Y3= NULL
        for (i in 1:n) {
            Y1[i] = rbinom(1, 1, Yp1[i])
            Y2[i] = rbinom(1, 1, Yp2[i])
            Y3[i] = rbinom(1, 1, Yp3[i])
        }
        Y = cbind(Y1, Y2, Y3)
        YW = cbind(Y, W)

        #observed outcomes
        Yobs = apply(YW, 1, function(x) x[1:3][x[4]])

        #table(Yobs, W)[2,]/table(W)*100
    }
    
    #==================================#
    #========= Extremely strong overlap===========#
    #==================================#
    if (overlap == "extremely strong") {
      Xcon = matrix(NA, nrow=n, ncol=5)
      Xcat = matrix(NA, nrow=n, ncol=5)
      
      for (i in 1:n) {
        Xcon[i,] = matrix(rnorm(5, mean = 0, sd = 1-0.01*W[i]),ncol=5);
        Xcat[i,] = matrix(sample (c(1,2,3), size = 5,replace =T, prob = c(.3,.3,.4)),ncol=5)
      }
      x1 = Xcon[,1]; x2 = Xcon[,2]; x3 = Xcon[,3]; x4 = Xcon[,4]; x5 = Xcon[,5];
      x6 = Xcat[,1]; x7 = Xcat[,2]; x8 = Xcat[,3]; x9 = Xcat[,4]; x10 = Xcat[,5]
      trtdat<- data.frame(W, Xcon, Xcat)
      
      #parallel response surface model
      #set.seed(3242019)
      tau1 = -4.7; tau2 = -4.5; tau3 = -5.0
      Yp1 = expit(tau1 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                  +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)
      Yp2 = expit(tau2 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                  +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)
      Yp3 = expit(tau3 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                  +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)
      
      #potential outcomes
      Y1 = Y2 = Y3= NULL
      for (i in 1:n) {
        Y1[i] = rbinom(1, 1, Yp1[i])
        Y2[i] = rbinom(1, 1, Yp2[i])
        Y3[i] = rbinom(1, 1, Yp3[i])
      }
      Y = cbind(Y1, Y2, Y3)
      YW = cbind(Y, W)
      
      #observed outcomes
      Yobs = apply(YW, 1, function(x) x[1:3][x[4]])
      
      #table(Yobs, W)[2,]/table(W)*100
    }

    
    #==================================#
    #=========strong overlap===========#
    #==================================#
    if (overlap == "strong") {
        Xcon = matrix(NA, nrow=n, ncol=5)
        Xcat = matrix(NA, nrow=n, ncol=5)

        for (i in 1:n) {
          Xcon[i,] = matrix(rnorm(5, mean = 0.05*W[i], sd = 1- 0.05*W[i]),ncol=5);
          Xcat[i,] = matrix(sample (c(1,2,3), size = 5,replace =T, prob = c(.3-0.001*W[i],.3+.001*W[i],.4)),ncol=5)
        }
        x1 = Xcon[,1]; x2 = Xcon[,2]; x3 = Xcon[,3]; x4 = Xcon[,4]; x5 = Xcon[,5];
        x6 = Xcat[,1]; x7 = Xcat[,2]; x8 = Xcat[,3]; x9 = Xcat[,4]; x10 = Xcat[,5]
        trtdat<- data.frame(W, Xcon, Xcat)

        #parallel response surface model
        #set.seed(3242019)
        tau1 = -4.7; tau2 = -4.5; tau3 = -5.0
        Yp1 = expit(tau1 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)
        Yp2 = expit(tau2 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)
        Yp3 = expit(tau3 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)

        #potential outcomes
        Y1 = Y2 = Y3= NULL
        for (i in 1:n) {
            Y1[i] = rbinom(1, 1, Yp1[i])
            Y2[i] = rbinom(1, 1, Yp2[i])
            Y3[i] = rbinom(1, 1, Yp3[i])
        }
        Y = cbind(Y1, Y2, Y3)
        YW = cbind(Y, W)

        #observed outcomes
        Yobs = apply(YW, 1, function(x) x[1:3][x[4]])

        #table(Yobs, W)[2,]/table(W)*100
    }

    #==================================#
    #========moderate overlap==========#
    #==================================#
    if (overlap == "moderate") {
        Xcon = matrix(NA, nrow=n, ncol=5)
        Xcat = matrix(NA, nrow=n, ncol=5)
        for (i in 1:n) {
          Xcon[i,] = matrix(rnorm(5,0.5,1)*(W[i]==1) + rnorm(5, 1, 1)*(W[i]==2) + rnorm(5, 1.5 , 1)*(W[i]==3),ncol=5);
          Xcat[i,] = matrix(sample (c(1,2,3), size = 5,replace =T, prob = c(.3-0.001*W[i],.3+.001*W[i],.4)),ncol=5)
        }
        x1 = Xcon[,1]; x2 = Xcon[,2]; x3 = Xcon[,3]; x4 = Xcon[,4]; x5  = Xcon[,5];
        x6 = Xcat[,1]; x7 = Xcat[,2]; x8 = Xcat[,3]; x9 = Xcat[,4]; x10 = Xcat[,5]
        trtdat<- data.frame(W, Xcon, Xcat)


        if (all_confounder == T) {
            tau1 = -6.0; tau2 = -4.5; tau3 = -5.0
            Yp1 = expit(tau1 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)
            Yp2 = expit(tau2 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)
            Yp3 = expit(tau3 + 0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
                +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 - .1*x10^2)

            #potential outcomes
            Y1 = Y2 = Y3= NULL
            for (i in 1:n) {
                Y1[i] = rbinom(1, 1, Yp1[i])
                Y2[i] = rbinom(1, 1, Yp2[i])
                Y3[i] = rbinom(1, 1, Yp3[i])
            }
            # observed outcomes
            Y = cbind(Y1,Y2,Y3)
            YW = cbind(Y,W)
            Yobs = apply(YW, 1, function(x) x[1:3][x[4]])

            #table(Yobs, W)[2,]/table(W)*100
        }

        if (all_confounder == F) {#only x1 and x5 are true confounnders
            tau1 = -3.0; tau2 = -4.8; tau3 = -6.2
            Yp1 = expit(tau1 + 0.2*x1 + 0.3*x6 - 0.2*x7 - 0.5*x8 + 0.6*x5 +.5*x1^2 - 0.1*x6^2 - 0.5*x7^2 + .15*x7^3 -  0.5*x1*x5)
            Yp2 = expit(tau2 + 0.2*x1 + 0.3*x6 - 0.2*x7 - 0.5*x8 + 0.6*x5 +.5*x1^2 - 0.1*x6^2 - 0.5*x7^2 + .15*x7^3 -  0.5*x1*x5)
            Yp3 = expit(tau3 + 0.2*x1 + 0.3*x6 - 0.2*x7 - 0.5*x8 + 0.6*x5 +.5*x1^2 - 0.1*x6^2 - 0.5*x7^2 + .15*x7^3 -  0.5*x1*x5)

            #potential outcomes
            Y1 = Y2 = Y3= NULL
            for (i in 1:n) {
                Y1[i] = rbinom(1, 1, Yp1[i])
                Y2[i] = rbinom(1, 1, Yp2[i])
                Y3[i] = rbinom(1, 1, Yp3[i])
            }

            # observed outcomes
            Y = cbind(Y1,Y2,Y3)
            YW = cbind(Y,W)
            Yobs = apply(YW, 1, function(x) x[1:3][x[4]])

            #table(Yobs, W)[2,]/table(W)*100
        }
    }

    #need the following two lines if run GBM; ignore if run BART
    colnames(trtdat) = c("trt_ind", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")
    trtdat$trt_ind = as.factor(trtdat$trt_ind)

    ###########################################################
    # true treatment effects
    # ATE(1,2), ATE(1,3), ATE(2,3)
    ATE12_RD = mean(Y[,1]) - mean(Y[,2])
    ATE13_RD = mean(Y[,1]) - mean(Y[,3])
    ATE23_RD = mean(Y[,2]) - mean(Y[,3])
    ATE12_RR = mean(Y[,1]) / mean(Y[,2])
    ATE13_RR = mean(Y[,1]) / mean(Y[,3])
    ATE23_RR = mean(Y[,2]) / mean(Y[,3])
    ATE12_OR = mean(Y[,1])/(1-mean(Y[,1])) / (mean(Y[,2])/(1-mean(Y[,2])))
    ATE13_OR = mean(Y[,1])/(1-mean(Y[,1])) / (mean(Y[,3])/(1-mean(Y[,3])))
    ATE23_OR = mean(Y[,2])/(1-mean(Y[,2])) / (mean(Y[,3])/(1-mean(Y[,3])))

    # ATT(1,2), ATT(1,3)
    ATT12_RD = mean(Y[W==1,1]) - mean(Y[W==1,2])
    ATT13_RD = mean(Y[W==1,1]) - mean(Y[W==1,3])
    ATT12_RR = mean(Y[W==1,1]) / mean(Y[W==1,2])
    ATT13_RR = mean(Y[W==1,1]) / mean(Y[W==1,3])
    ATT12_OR = mean(Y[W==1,1])/(1-mean(Y[W==1,1])) / (mean(Y[W==1,2])/(1-mean(Y[W==1,2])))
    ATT13_OR = mean(Y[W==1,1])/(1-mean(Y[W==1,1])) / (mean(Y[W==1,3])/(1-mean(Y[W==1,3])))

    Ests = list(ATE12_RD=ATE12_RD, ATE13_RD=ATE13_RD, ATE23_RD=ATE23_RD,
                ATE12_RR=ATE12_RR, ATE13_RR=ATE13_RR, ATE23_RR=ATE23_RR,
                ATE12_OR=ATE12_OR, ATE13_OR=ATE13_OR, ATE23_OR=ATE23_OR,
                ATT12_RD=ATT12_RD, ATT13_RD=ATT13_RD, ATT12_RR=ATT12_RR,
                ATT13_RR=ATT13_RR, ATT12_OR=ATT12_OR, ATT13_OR=ATT13_OR)

    return(list(trtdat=trtdat, Y=Y, Yobs=Yobs, n=n, Ests=Ests))
}


set.seed(3242019)
# mydata_p3_scen1 = data_gen(n = 4e5, overlap = "weak")
# save(mydata_p3_scen1, file="Rdata/sim_p3_scen1.RData")
# mydata_p3_scen2 = data_gen(n = 4e5, overlap = "strong")
# save(mydata_p3_scen2, file="Rdata/sim_p3_scen2.RData")
# mydata_p3_scen3 = data_gen(n = 4e5, overlap = "moderate")
# save(mydata_p3_scen3, file="Rdata/sim_p3_scen3.RData")
# mydata_p3_scen4 = data_gen(n = 4e5, overlap = "extremely strong")
# save(mydata_p3_scen4, file="Rdata/sim_p3_scen4.RData")
