#============================================================================#
# Simulation design - part 1                                                 #
# sample size 1500, 4000, 9500 with ratio of units 1:1:1, 1:4:3 and 1:10:8   #
# number of confounders 10                                                   #
# Author: Liangyuan Hu & Chenyang Gu                                         #
# Feb 11, 2020                                                            #
#============================================================================#

expit = function(x) {exp(x)/(1+exp(x))}

data_gen = function(n = 9500, ratio = 3) {

    p = 10

    #treatment assignment model
    if (ratio == 1) {alpha1 = 0.32;  alpha2 = 0.28}
    if (ratio == 2) {alpha1 = -0.9;  alpha2 = 0.7}
    if (ratio == 3) {alpha1 = -2.4; alpha2 = 0.52}

    #set.seed(3242019)
    #continuous covariates
    Xcon = matrix(rnorm(p*n), nrow=n, ncol=p/2)
    x1 = Xcon[,1]; x2 = Xcon[,2]; x3 = Xcon[,3]; x4 = Xcon[,4]; x5 = Xcon[,5]
    #categorical covariates
    Xcat = matrix(sample(0:2, n*p, replace = T, prob = c(.3,.3,.4)), nrow=n, ncol=p/2)
    x6 = Xcat[,1]; x7 = Xcat[,2]; x8 = Xcat[,3]; x9 = Xcat[,4]; x10 = Xcat[,5]

    ex1 = exp(alpha1 + .2*x1 + .4*x2 + .3*x3 + .4*x4 + .1*x5 + .2*x6 - .8*x7 - 1.1*x8 + .5*x9 + 0.5*x10
                + .2*x1^2 + .3*x2^2 + .4*x4^2 + .8*x1*x2 + .5*x1*x6 + .4*x4*x8)
    ex2 = exp(alpha2 + .5*x1 + .8*x2 + .6*x3 + .2*x4 + .25*x5 + .4*x6 - 1.2*x7 - 1.5*x8 - .3*x9 + 1.5*x10
                + .2*x1^2 + .7*x2^2 + .2*x4^2 + .25*x1*x2 + .3*x1*x6 + .6*x4*x8)

    Wp1 = ex1 / (1 + ex1 + ex2)
    Wp2 = ex2 / (1 + ex1 + ex2)
    Wp3 = 1 - Wp1 - Wp2

    W = NULL
    for (i in 1:n) W[i] = sample(c(1,2,3), size=1, replace=T, prob=c(Wp1[i],Wp2[i],Wp3[i]))
    #table(W)[2]/table(W)[1]
    #table(W)[3]/table(W)[1]

    trtdat = data.frame(W, Xcon, Xcat)
    colnames(trtdat) = c("trt_ind", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")
    trtdat$trt_ind = as.factor(trtdat$trt_ind)
    #require(nnet)
    #summary( multinom(W ~ ., data = trtdat))

    #parallel response surface model
    #set.seed(3242019)
    tau1 = -4.5; tau2 = -4.7; tau3 = -5.0
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
    #table(Y1)[2]/sum(table(Y1))
    #table(Y2)[2]/sum(table(Y2))
    #table(Y3)[2]/sum(table(Y3))
    Y = cbind(Y1, Y2, Y3)
    YW = cbind(Y, W)
 
    #observed outcomes
    Yobs = apply(YW, 1, function(x) x[1:3][x[4]])
    #table(Yobs)[2]/(sum(table(Yobs)))

   
    #table(Yobs, W)[2,]/table(W)*100

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


# set.seed(3242019)
# mydata_p1_scen1 = data_gen(n = 4e5, ratio =1)
# save(mydata_p1_scen1, file="Rdata/sim_p1_scen1.RData")
# mydata_p1_scen2 = data_gen(n = 4e5, ratio =2)
# save(mydata_p1_scen2, file="Rdata/sim_p1_scen2.RData")
# mydata_p1_scen3 = data_gen(n = 4e5, ratio = 3)
# save(mydata_p1_scen3, file="Rdata/sim_p1_scen3.RData")
