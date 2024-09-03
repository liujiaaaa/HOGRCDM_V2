#library(gtools)
INI_CC_Lognormal<-function(DATA, N,J,K,eplison_Layer0=1e-04){
  
  
  
  Init_B_Intec1<-colMeans(DATA)
  DATA1<-scale(DATA,center=T,scale=F)
  
  
  Svd_Layer0_SEC<-svd(DATA1)
  singular_Layer0_SEC_left<-Svd_Layer0_SEC$u
  singular_Layer0_SEC_right<-Svd_Layer0_SEC$v
  singular_Layer0_SEC_middle<-Svd_Layer0_SEC$d
  
  
  
  #  K_tuta<-max(K+1,max(which(singular_Layer0_SEC_middle>=1.01*sqrt(N))))
  # Slope<-(1/sqrt(N))*singular_Layer0_SEC_right[,1:K_tuta]%*%diag(singular_Layer0_SEC_middle[1:K_tuta])
  # 
  # Init_B_Intec1<-Slope[,1]
  # theta<-sqrt(N)*singular_Layer0_SEC_left[,1:K]
  # Slope1<-varimax(Slope[,-1], normalize = T, eps = 1e-5)
  # Init_B_Slope1<-as.matrix((Slope1$loadings))
  # Init_B_Slope1 <- as.matrix(Matrix::Matrix(Init_B_Slope1))
  
  
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
  
  
  B1<-t(solve(t(A_hat)%*%A_hat)%*%t(A_hat)%*%(Z0_hat))*G1
  for(tt in 1: ncol(B1)){
    if(B1[,tt][which.max(abs(B1[,tt]))]<0) B1[,tt]<--B1[,tt]
  }
  
  return(list(Init_B_Slope1=B1,Init_B_Intec1=Init_B_Intec1, A_hat=A_hat, G1= G1))
}
