GenerateData<-function(SizeList,ModelSetList,TrueP_List_Bott,TrueP_List_High,Q_B){
  
  
  N=SizeList$N
  J=SizeList$J
  K=SizeList$K
  D=SizeList$D
  L=2^K
  
  Q_Aug=TrueP_List_Bott$Q_Aug
 
  
  
  
  
  if(!is.null(TrueP_List_High$thetaT)){
    thetaT<-TrueP_List_High$thetaT
  }else{
    thetaT<-mvrnorm(n=N,TrueP_List_High$Mu_thetaT,TrueP_List_High$Sigma_thetaT)
  }
  
  
  
  
  
  AlphaT<-generate4PNO11(A=TrueP_List_High$Slope_H,b= TrueP_List_High$Interc_H,
                         c=rep(0,K),d=rep(1,K),theta=thetaT,type="4MPNO",
                         settings=SizeList)$U
  
  
  
  if(ModelSetList$ModelStruc_Bottom=="Main_effect"){
    K_q<-K
    AlphaT_Aug<- AlphaT
    
  }else if(ModelSetList$ModelStruc_Bottom=="All_effect"){
    K_q<-L
    AlphaT_Aug<- AugmenInteraction(OrigriData=AlphaT,K)$OrigriData_Aug
  }else if(ModelSetList$ModelStruc_Bottom=="DINA"){
    K_q<-L
    AlphaT_Aug<- AugmenInteraction(OrigriData=AlphaT,K)$OrigriData_Aug
  }
  
  
  if(ModelSetList$ModelDistri_Bottom=="Lognormal"){
    Alpha_J<-list()
    Pro_Response<-matrix(NA,N,J)
    for(j in 1:J){
      Alpha_J[[j]]<-AlphaT_Aug[,Q_Aug[j,]!=0]
      Pro_Response[,j]<-rep(TrueP_List_Bott$Interc_B[j],N)+
        rowSums((rep(1,N)%*%t(TrueP_List_Bott$Slope_B[j,][Q_Aug[j,]!=0]))*as.matrix(Alpha_J[[j]]))
    }
    
    sd_M<-rep(1,N)%*%t(TrueP_List_Bott$sd)
    Res1<-rnorm(N*J,mean=as.vector(Pro_Response),sd=sd_M)
    Res<-matrix(exp(Res1),N,J)
  }
  
  
  ######## Complete other model setting lately
  
 
  return(Res=Res)
  
}