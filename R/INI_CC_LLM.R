#library(gtools)
INI_CC_LLM<-function(DATA, N,J,K,eplison_Layer0=1e-04){
  Res<-DATA
  K1<-K
  
  Svd_Layer0<-svd(Res)
  singular_Layer0_left<-Svd_Layer0$u
  singular_Layer0_right<-Svd_Layer0$v
  singular_Layer0_middle<-Svd_Layer0$d
  Te<-which(singular_Layer0_middle>=1.01*sqrt(N))
  
  if(length(Te)==0){
    MAX1<-which.max(singular_Layer0_middle)
  }else{
    MAX1<-max(Te)
  }
  if(is.infinite((  MAX1)))   MAX1<-0
  K_tuta<-max(K1+1,MAX1) 
  X_Layer0<-singular_Layer0_left[,1:K_tuta]%*%diag(singular_Layer0_middle[1:K_tuta])%*%
    t(singular_Layer0_right[,1:K_tuta])
  
  
  X_Layer0_hat<-matrix(NA,N,J)
  X_Layer0_hat[(X_Layer0>eplison_Layer0&X_Layer0<1-eplison_Layer0)]<-
    X_Layer0[(X_Layer0>eplison_Layer0&X_Layer0<1-eplison_Layer0)]
  
  
  X_Layer0_hat[X_Layer0<eplison_Layer0]<-eplison_Layer0
  X_Layer0_hat[X_Layer0>1-eplison_Layer0]<-1-eplison_Layer0
  ######################################
  M_tuta_Layer0<-gtools::logit(X_Layer0_hat)
  #########################################
  Init_B_Intec1<-colMeans(M_tuta_Layer0)
  
  
  M_hat_Layer0<-scale(M_tuta_Layer0,center=T,scale=F)
  
  
  
  Svd_Layer0_SEC<-svd(M_hat_Layer0)
  singular_Layer0_SEC_left<-Svd_Layer0_SEC$u
  singular_Layer0_SEC_right<-Svd_Layer0_SEC$v
  singular_Layer0_SEC_middle<-Svd_Layer0_SEC$d
  
  
  
  Z0_hat<- singular_Layer0_SEC_left[,1:K]%*%diag(singular_Layer0_SEC_middle[1:K])%*%t( singular_Layer0_SEC_right[,1:K])
  
  V_tuta<-varimax(singular_Layer0_SEC_right[,1:K], normalize =F, eps = 1e-5)$loadings
  
  # V_tuta<- targetQ(singular_Layer0_SEC_right[,1:K], Target=TARGET)$loadings
  V_tuta[abs(V_tuta)<1/(2*sqrt(J))]<-0
  
  G1<- as.matrix(Matrix::Matrix( V_tuta))
  G1[ G1!=0]<-1
  
  
  
  A_0<-Z0_hat%*% V_tuta%*%solve(t(V_tuta)%*% V_tuta)
  A_hat<- A_0
  A_hat[ A_hat>0]<-1
  A_hat[ A_hat<0| A_hat==0]<-0
  
  
  B1<-1/2*t(solve(t(A_hat)%*%A_hat)%*%t(A_hat)%*%(Z0_hat))*G1
  for(tt in 1: ncol(B1)){
    if(B1[,tt][which.max(abs(B1[,tt]))]<0) B1[,tt]<--B1[,tt]
  }
  
  
  
  return(list(Init_B_Slope1=B1,Init_B_Intec1= Init_B_Intec1, A_hat= A_hat,G1=G1))
}
