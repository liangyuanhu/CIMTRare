#######################################################################
# Functions to estimate causal effects of multiple treatment using BART
# Author: Chenyang Gu
# Date: 02/23/2019
#######################################################################


# Estimate average treatment effect (ATE)
bart.multiTrt.ate = function(y, x, trt, k=2, ntree=100, ndpost=1000, nskip=1000) {
    
    n1 = sum(trt==1)
    n2 = sum(trt==2)
    n3 = sum(trt==3)
    
    xt = cbind(trt,x)
    
    
    # Fit BART
    bart_mod = pbart(x.train = xt, y.train = y, k = k, ntree = ntree, ndpost = ndpost, nskip = nskip)
    
    
    # Predict potential outcomes for trt=1
    xp1 = xt[trt==1,]
    xp2 = xp1
    xp3 = xp1
    xp2[,1] = 2  # switch treatment label 1 to 2
    xp3[,1] = 3  # switch treatment label 1 to 3
    
    bart_pred11 = pwbart(xp1, bart_mod$treedraws)
    bart_pred12 = pwbart(xp2, bart_mod$treedraws)
    bart_pred13 = pwbart(xp3, bart_mod$treedraws)
    
    pred_prop11 = pnorm(bart_pred11)
    pred_prop12 = pnorm(bart_pred12)
    pred_prop13 = pnorm(bart_pred13)
    
    
    # Predict potential outcomes for trt=2
    xp2 = xt[trt==2,]
    xp1 = xp2
    xp3 = xp2
    xp1[,1] = 1  # switch treatment label 2 to 1
    xp3[,1] = 3  # switch treatment label 2 to 3
    
    bart_pred21 = pwbart(xp1, bart_mod$treedraws)
    bart_pred22 = pwbart(xp2, bart_mod$treedraws)
    bart_pred23 = pwbart(xp3, bart_mod$treedraws)
    
    pred_prop21 = pnorm(bart_pred21)
    pred_prop22 = pnorm(bart_pred22)
    pred_prop23 = pnorm(bart_pred23)
    
    # Predict potential outcomes for trt=3
    xp3 = xt[trt==3,]
    xp1 = xp3
    xp2 = xp3
    xp1[,1] = 1  # switch treatment label 3 to 1
    xp2[,1] = 2  # switch treatment label 3 to 2
    
    bart_pred31 = pwbart(xp1, bart_mod$treedraws)
    bart_pred32 = pwbart(xp2, bart_mod$treedraws)
    bart_pred33 = pwbart(xp3, bart_mod$treedraws)
    
    pred_prop31 = pnorm(bart_pred31)
    pred_prop32 = pnorm(bart_pred32)
    pred_prop33 = pnorm(bart_pred33)
    
    
    # Estimate causal effects
    RD12.est = RR12.est = OR12.est = NULL
    RD13.est = RR13.est = OR13.est = NULL
    RD23.est = RR23.est = OR23.est = NULL
    
    for (m in 1:ndpost) {
        
        # Estimate E(Y1), E(Y2), E(Y3)
        y1 = c(rbinom(n1, 1, pred_prop11[m,]), rbinom(n2, 1, pred_prop21[m,]), rbinom(n3, 1, pred_prop31[m,]))
        y2 = c(rbinom(n1, 1, pred_prop12[m,]), rbinom(n2, 1, pred_prop22[m,]), rbinom(n3, 1, pred_prop32[m,]))
        y3 = c(rbinom(n1, 1, pred_prop13[m,]), rbinom(n2, 1, pred_prop23[m,]), rbinom(n3, 1, pred_prop33[m,]))
        
        y1.pred = mean(y1)
        y2.pred = mean(y2)
        y3.pred = mean(y3)
        
        # Calculate risk difference (RD)
        RD12.est[m] = y1.pred - y2.pred
        RD13.est[m] = y1.pred - y3.pred
        RD23.est[m] = y2.pred - y3.pred
        
        # Calculate relative risk (RR)
        RR12.est[m] = y1.pred / y2.pred
        RR13.est[m] = y1.pred / y3.pred
        RR23.est[m] = y2.pred / y3.pred
        
        # Calculate  odds ratio (OR)
        OR12.est[m] = (y1.pred / (1 - y1.pred)) / (y2.pred / (1 - y2.pred))
        OR13.est[m] = (y1.pred / (1 - y1.pred)) / (y3.pred / (1 - y3.pred))
        OR23.est[m] = (y2.pred / (1 - y2.pred)) / (y3.pred / (1 - y3.pred))
    }
    
    ate12 = postSumm(RD12.est, RR12.est, OR12.est)
    ate13 = postSumm(RD13.est, RR13.est, OR13.est)
    ate23 = postSumm(RD23.est, RR23.est, OR23.est)
    
    list(ATE12 = round(ate12, digits=3),
         ATE13 = round(ate13, digits=3),
         ATE23 = round(ate23, digits=3))
}



# Estimate average treatment effect on the treated (ATT)
# The default reference group is 1st group
bart.multiTrt.att = function(y, x, trt, k=2, ntree=100, ndpost=1000, nskip=1000) {
    
    n1 = sum(trt==1)
    n2 = sum(trt==2)
    n3 = sum(trt==3)
    
    xt = cbind(trt,x)
    
    # Fit BART
    bart_mod = pbart(x.train = xt, y.train = y, k = k, ntree = ntree, ndpost = ndpost, nskip = nskip)
    
    
    # Predict potential outcomes for trt=1
    xp1 = xt[trt==1,]
    xp2 = xp1
    xp3 = xp1
    xp2[,1] = 2  # switch treatment label 1 to 2
    xp3[,1] = 3  # switch treatment label 1 to 3
   
    bart_pred11 = pwbart(xp1, bart_mod$treedraws)
    bart_pred12 = pwbart(xp2, bart_mod$treedraws)
    bart_pred13 = pwbart(xp3, bart_mod$treedraws)
    
    pred_prop11 = pnorm(bart_pred11)
    pred_prop12 = pnorm(bart_pred12)
    pred_prop13 = pnorm(bart_pred13)
    
    
    # Estimate causal effects
    RD12.est = RR12.est = OR12.est = NULL
    RD13.est = RR13.est = OR13.est = NULL
    
    for (m in 1:ndpost) {
        
        # Estimate E(Y1|trt=1), E(Y2|trt=1), E(Y3|trt=1)
        y1.pred = mean(rbinom(n1, 1, pred_prop11))
        y2.pred = mean(rbinom(n1, 1, pred_prop12))
        y3.pred = mean(rbinom(n1, 1, pred_prop13))
        
        # Calculate risk difference (RD)
        RD12.est[m] = y1.pred - y2.pred
        RD13.est[m] = y1.pred - y3.pred
        
        # Calculate relative risk (RR)
        RR12.est[m] = y1.pred / y2.pred
        RR13.est[m] = y1.pred / y3.pred
        
        # Calculate odds ratio (OR)
        OR12.est[m] = (y1.pred / (1 - y1.pred)) / (y2.pred / (1 - y2.pred))
        OR13.est[m] = (y1.pred / (1 - y1.pred)) / (y3.pred / (1 - y3.pred))
    }
    
    att12 = postSumm(RD12.est, RR12.est, OR12.est)
    att13 = postSumm(RD13.est, RR13.est, OR13.est)
    
    list(ATT12 = round(att12, digits=3),
         ATT13 = round(att13, digits=3))
}


# Use BART to estimate causal effects of multiple treatments
bart.multiTrt = function(y, x, trt, estimand="ATE", k=2, ntree=100, ndpost=1000, nskip=1000) {
    
    # Data structure
    #        Y(1) Y(2) Y(3)
    # trt=1   *    ?    ?
    # trt=2   ?    *    ?
    # trt=3   ?    ?    *
    
    #        Y(1) Y(2) Y(3)
    # trt=1  y11  y12  y13
    # trt=2  y21  y22  y23
    # trt=3  y31  y32  y33
    
    if (estimand=="ATE") {
        bart.est = bart.multiTrt.ate(y, x, trt, k=2, ntree=100, ndpost=1000, nskip=1000)
    }
    
    if (estimand=="ATT") {
        bart.est = bart.multiTrt.att(y, x, trt, k=2, ntree=100, ndpost=1000, nskip=1000)
    }
    
    return(bart.est)
}



























