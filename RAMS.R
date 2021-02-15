#============================================================================#
# RAMS                                                                       #
# Author: Chenyang Gu                                                        #                                                               
#============================================================================#

RegSpline = function(y, trt, ps1, ps2) {
  
  n = length(trt)
  
  # logit of propensity scores
  logit_ps1 = qlogis(ps1)
  logit_ps2 = qlogis(ps2)
  mod.spline = gam(y ~ trt + te(logit_ps1,logit_ps2), family = binomial(link="logit"))
  
  # prediction
  # predict potential outcomes Y(1)
  newdata1 = data.frame(trt=rep(1,n), logit_ps1=logit_ps1, logit_ps2=logit_ps2)
  spline.pred1 = plogis(predict(mod.spline, newdata = newdata1))
  
  # predict potential outcomes Y(2)
  newdata2 = data.frame(trt=rep(2,n), logit_ps1=logit_ps1, logit_ps2=logit_ps2)
  spline.pred2 = plogis(predict(mod.spline, newdata = newdata2))
  
  # predict potential outcomes Y(3)
  newdata3 = data.frame(trt=rep(3,n), logit_ps1=logit_ps1, logit_ps2=logit_ps2)
  spline.pred3 = plogis(predict(mod.spline, newdata = newdata3))
  
  y1.hat = mean(spline.pred1)
  y2.hat = mean(spline.pred2)
  y3.hat = mean(spline.pred3)
  
  # # relative risk
  # ate12 = y1.hat / y2.hat
  # ate13 = y1.hat / y3.hat
  # ate23 = y2.hat / y3.hat
  
  # relative difference
  ate12 = y1.hat - y2.hat
  ate13 = y1.hat - y3.hat
  ate23 = y2.hat - y3.hat
  
  res = c(ate12, ate13, ate23)
  return(res)
}