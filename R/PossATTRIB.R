PossATTRIB<-function(K){
  L<-2^K
  Poss_Attr<-matrix(NA,2^K,K)
  for(l in 1:K){
    repTimes<-2^(l-1)
    Poss_Attr[,l]<-rep(c(rep(1,L/2^l),rep(0,L/2^l)),repTimes)
  }
  Poss_Attr
  
}