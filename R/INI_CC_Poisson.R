#library(gtools)
INI_CC_Poisson<-function(DATA, N,J,K,scale=1){
  Res<-log(DATA+1) #Res>= log(1)
  K1<-K
  eplison_Layer0<-log(1)
  
  Svd_Layer0<-svd(Res)
  singular_Layer0_left<-Svd_Layer0$u
  singular_Layer0_right<-Svd_Layer0$v
  singular_Layer0_middle<-Svd_Layer0$d
  MAX1<-max(which(singular_Layer0_middle>=1.01*sqrt(N)))
  if(is.infinite((  MAX1)))   MAX1<-0
  K_tuta<-max(K1+1,MAX1)
  X_Layer0<-singular_Layer0_left[,1:K_tuta]%*%diag(singular_Layer0_middle[1:K_tuta])%*%
    t(singular_Layer0_right[,1:K_tuta])
  
  
  X_Layer0_hat<-matrix(NA,N,J)
  
  X_Layer0_hat[(X_Layer0>eplison_Layer0)]<-
    X_Layer0[(X_Layer0>eplison_Layer0)]
  
  X_Layer0_hat[(X_Layer0==eplison_Layer0)]<-
    X_Layer0[(X_Layer0==eplison_Layer0)]
  
  X_Layer0_hat[X_Layer0<eplison_Layer0]<-eplison_Layer0
  
  M_tuta_Layer0<-exp(X_Layer0_hat)-1
  
  Init_B_Intec1<-colMeans(as.matrix(M_tuta_Layer0))
  M_hat_Layer0<-scale(M_tuta_Layer0,center=T,scale=F)
  
  
  
  Svd_Layer0_SEC<-svd(M_hat_Layer0)
  singular_Layer0_SEC_left<-Svd_Layer0_SEC$u
  singular_Layer0_SEC_right<-Svd_Layer0_SEC$v
  singular_Layer0_SEC_middle<-Svd_Layer0_SEC$d
  
  
  
  
  
  V_tuta<-varimax(singular_Layer0_SEC_right[,1:K], normalize =F, eps = 1e-5)$loadings
  
  # V_tuta<- targetQ(singular_Layer0_SEC_right[,1:K], Target=TARGET)$loadings
  V_tuta[abs(V_tuta)<1/(2*sqrt(J))]<-0
  
  G1<- as.matrix(Matrix::Matrix( V_tuta))
  G1[ G1!=0]<-1
  
  
  
  A_0<-M_hat_Layer0%*% V_tuta%*%solve(t(V_tuta)%*% V_tuta)
  A_hat<- A_0
  A_hat[ A_hat>0]<-1
  A_hat[ A_hat<0| A_hat==0]<-0
  
  
  B1<-t(solve(t(A_hat)%*%A_hat)%*%t(A_hat)%*%(M_hat_Layer0))*G1/scale
  for(tt in 1: ncol(B1)){
    if(B1[,tt][which.max(abs(B1[,tt]))]<0) B1[,tt]<--B1[,tt]
  }
  
  Init_B_Intec1<-Init_B_Intec1/scale
  return(list(Init_B_Slope1=B1,Init_B_Intec1=Init_B_Intec1,A_hat=A_hat, G1=G1,V_tuta=V_tuta))
}
