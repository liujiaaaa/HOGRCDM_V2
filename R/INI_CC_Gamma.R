#library(gtools)
INI_CC_Gamma<-function(DATA, N,J,K,scale=1){

  
  
  mm<-apply(Res,2,mean)
  
  ss<-apply(Res,2,var)
  
  Shape<-mm^2/ss
  
  Rate<-mm/ss
  
  
  

  # 计算 Gamma 分布的均值和方差
  mu_gamma <-apply(Res,2,mean)
  sigma2_gamma <- apply(Res,2,var)
  
  # 计算对数正态分布的对数均值和对数标准差
  mu_log <- log(mu_gamma) - 0.5 * log(1 + sigma2_gamma / mu_gamma^2)
  sigma_log <- sqrt(log(1 + sigma2_gamma / mu_gamma^2))
  
  # 假设你有一组服从 Gamma 分布的数据 X_Gamma
  X_Gamma <- Res
  
  # 将 Gamma 分布数据转换为标准化的数据
  Z_Gamma <- apply( X_Gamma,2, scale)
  
  # 
  # X_LogNormal <- log(mu_log + Z_Gamma * sigma_log)
  
  
  X_Normal <- (mu_log + Z_Gamma * sigma_log)
  
 
  
  IZ<-INI_CC_Lognormal(DATA= X_Normal, N,J,K,eplison_Layer0=1e-04)
  
  Init_B_Slope1<-IZ$Init_B_Slope1
  Init_B_Intec1<-IZ$Init_B_Intec1
  A_hat<-IZ$A_hat
  G1<-IZ$G1

  
  
  return(list(Init_B_Slope1=Init_B_Slope1,Init_B_Intec1=Init_B_Intec1,A_hat=A_hat, G1=G1,
              Shape=Shape,Rate=Rate))
}
