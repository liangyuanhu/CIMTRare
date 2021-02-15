#============================================================================#
# Function to calculate ATE for binary outcome with multiple (3) treatments  #
# Author: Liangyuan Hu                                                       #                                                               
#============================================================================#

ate_fun<- function(wt1 = ate.wt.1, wt2 = ate.wt.2, wt3 = ate.wt.3){
  mu_1_hat.iptw = sum(y[trt_ind==1] * wt1[trt_ind==1]) / sum(wt1[trt_ind==1])
  mu_2_hat.iptw = sum(y[trt_ind==2] * wt2[trt_ind==2]) / sum(wt2[trt_ind==2])
  mu_3_hat.iptw = sum(y[trt_ind==3] * wt3[trt_ind==3]) / sum(wt3[trt_ind==3])
  RD12 = mu_1_hat.iptw - mu_2_hat.iptw            
  RD13 = mu_1_hat.iptw - mu_3_hat.iptw  
  RD23 = mu_2_hat.iptw - mu_3_hat.iptw
  RR12 = mu_1_hat.iptw/mu_2_hat.iptw
  RR13 = mu_1_hat.iptw/mu_3_hat.iptw
  RR23 = mu_2_hat.iptw/mu_3_hat.iptw
  OR12 = (mu_1_hat.iptw/(1-mu_1_hat.iptw))/(mu_2_hat.iptw/(1-mu_2_hat.iptw))
  OR13 = (mu_1_hat.iptw/(1-mu_1_hat.iptw))/(mu_3_hat.iptw/(1-mu_3_hat.iptw))
  OR23 = (mu_2_hat.iptw/(1-mu_2_hat.iptw))/(mu_3_hat.iptw/(1-mu_3_hat.iptw))
  res = list(RD12, RD13, RD23, RR12,RR13, RR23, OR12,OR13, OR23)
  return (res)
}
trunc_fun<-function(x){pmin(quantile(x,.95),pmax(quantile(x,.05),x))}