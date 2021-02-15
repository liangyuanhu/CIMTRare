#######################################################################
# Function to summarize posterior samples of RD, RR and OR
# Author: Chenyang Gu
# Date: 02/23/2019
#######################################################################


postSumm = function(RD.est, RR.est, OR.est) {
    # Risk difference (RD)
    RD.mean = mean(RD.est)
    RD.se = sd(RD.est)
    RD.lower = quantile(RD.est, probs=0.025, na.rm = T)
    RD.upper = quantile(RD.est, probs=0.975, na.rm = T)
    
    # Relative risk (RR)
    RR.mean = mean(RR.est)
    RR.se = sd(RR.est)
    RR.lower = quantile(RR.est, probs=0.025, na.rm = T)
    RR.upper = quantile(RR.est, probs=0.975, na.rm = T)
    
    # Odds ratio (OR)
    OR.mean = mean(OR.est)
    OR.se = sd(OR.est)
    OR.lower = quantile(OR.est, probs=0.025, na.rm = T)
    OR.upper = quantile(OR.est, probs=0.975, na.rm = T)
    
    # summarize results
    RD = c(RD.mean, RD.se, RD.lower, RD.upper)
    RR = c(RR.mean, RR.se, RR.lower, RR.upper)
    OR = c(OR.mean, OR.se, OR.lower, OR.upper)
    
    res = rbind(RD, RR, OR)
    colnames(res) = c("EST","SE","LOWER","UPPER")
    return(res)
}























