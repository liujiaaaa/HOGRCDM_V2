#Generate data with given parameters
generate4PNO11<-function(A,b,c,d=NULL,theta,type="4MPNO",settings){
 if(type=="4MPNO"&is.null(d)) stop("Parameter 'd' must be input for 4MPNO model")
  N<-settings$N
  J<-settings$K
  D<-settings$D
  if(D>1){
    TAb<-theta%*%A+rep(1,N)%*%t(b)
  }else{
    TAb<-theta%*%t(A)+rep(1,N)%*%t(b)
  }
 
  PHI1<-pnorm(TAb)
  #PHI2<-1-PHI1
  #Prob<-(rep(1,N)%*%t(d))*PHI1+(rep(1,N)%*%t(c))*PHI2
  if(type=="4MPNO"){
    Prob<-(rep(1,N)%*%t(c))+(rep(1,N)%*%t(d-c))*PHI1
  }else if(type=="3MPNO"){
    d<-NULL
    Prob<-(rep(1,N)%*%t(c))+(rep(1,N)%*%t(1-c))*PHI1
  } 
  ru<-matrix(runif(N*J,0,1),N,J)
  U<-matrix(as.numeric((Prob-ru)>=0),N,J)
  if(D>1){
    XI=cbind(t(A),b,c,d)
  }else{
    XI=cbind(A,b,c,d)
  }
 
  return(list(U=U,A=A,b=b,c=c,d=d,theta=theta,XI=XI))
}
