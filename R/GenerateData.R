GenerateData<-function(SizeList,ModelSetList,TrueP_List_Bott,TrueP_List_High,Q_B){
  
  
  N=SizeList$N
  J=SizeList$J
  K=SizeList$K
  D=SizeList$D
  L=2^K
  ModelStruc_Bottom<-ModelSetList$ModelStruc_Bottom
  
 
  
  
  
  
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
    
    Q_Aug<-Q_B
    
  }else if(ModelSetList$ModelStruc_Bottom=="All_effect"){
    K_q<-L
    AlphaT_Aug<- AugmenInteraction(OrigriData=AlphaT,K)$OrigriData_Aug
    
    Q_Aug<-Q_B
    q_num1<-K
    for(kk in 2:K){
      combinations <- t(combn(1:K, kk))
      te<-apply( Q_B,1,FUN= comFUN, combinations)
      # print(dim(as.matrix(te)))
      # print(choose(K1,kk))
      if(kk<K) te<-t(te)
      Q_Aug<-cbind(Q_Aug,te)
      q_num1<-q_num1+choose(K,kk)
    }
    
    
  }else if(ModelSetList$ModelStruc_Bottom=="DINA"){
    K_q<-L
    AlphaT_Aug<- AugmenInteraction(OrigriData=AlphaT,K)$OrigriData_Aug
    
    
    DINA_INDEX<-rep(1,J)
    Q_Aug<-Q_B
    q_num1<-K
    for(kk in 2:K){
      combinations <- t(combn(1:K, kk))
      te<-apply( Q_B,1,FUN= comFUN, combinations)
      # print(dim(as.matrix(te)))
      # print(choose(K1,kk))
      if(kk<K){
        te<-t(te)
        DINA_INDEX<-  (1-as.numeric(rowSums(te)>0))*DINA_INDEX+kk*as.numeric(rowSums(te)>0)
        
        
        
      }else if(kk==K){
        DINA_INDEX<-  (1-as.numeric((te)>0))*DINA_INDEX+kk*as.numeric((te)>0)
      }
      indd<-which(DINA_INDEX==kk)
      if(length(indd)>0)  Q_Aug[indd,]<-0
      Q_Aug<-cbind(Q_Aug,te)
      q_num1<-q_num1+choose(K,kk)
    }
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
    
  }else if(ModelSetList$ModelDistri_Bottom=="LLM"){
    Alpha_J<-list()
    Pro_Response<-matrix(NA,N,J)
    for(j in 1:J){
      Alpha_J[[j]]<-AlphaT_Aug[,Q_Aug[j,]!=0]
      Pro_Response[,j]<-rep(TrueP_List_Bott$Interc_B[j],N)+
        rowSums((rep(1,N)%*%t(TrueP_List_Bott$Slope_B[j,][Q_Aug[j,]!=0]))*as.matrix(Alpha_J[[j]]))
    }
    
    PPRES<-plogis( Pro_Response)
    rU<-runif(N*J,0,1)
    Res<-matrix(as.numeric(rU<PPRES),N,J) #Res is the observed response Data
    
  } else if(ModelSetList$ModelDistri_Bottom=="Poisson"){
    Alpha_J<-list()
    Pro_Response<-matrix(NA,N,J)
    for(j in 1:J){
      Alpha_J[[j]]<-AlphaT_Aug[,Q_Aug[j,]!=0]
      Pro_Response[,j]<-rep(TrueP_List_Bott$Interc_B[j],N)+
        rowSums((rep(1,N)%*%t(TrueP_List_Bott$Slope_B[j,][Q_Aug[j,]!=0]))*as.matrix(Alpha_J[[j]]))
    }
    
    Res1<-rpois(N*J,lambda=as.vector(Pro_Response))
    Res<-matrix(Res1,N,J)

  }else if (ModelSetList$ModelDistri_Bottom=="Gamma"){
    Alpha_J<-list()
    Pro_Response<-matrix(NA,N,J)
    for(j in 1:J){
      Alpha_J[[j]]<-AlphaT_Aug[,Q_Aug[j,]!=0]
      Pro_Response[,j]<-rep(TrueP_List_Bott$Interc_B[j],N)+
        rowSums((rep(1,N)%*%t(TrueP_List_Bott$Slope_B[j,][Q_Aug[j,]!=0]))*as.matrix(Alpha_J[[j]]))
    }
    ShapeM<-rep(1,N)%*%t(TrueP_List_Bott$Shape)
    Res1<- rgamma(N*J, shape=  ShapeM, rate = as.vector(Pro_Response))
    Res<-matrix( Res1,N,J)
    
  }
  
  
  ######## Complete other model setting lately
  
 
  return(Res=Res)
  
}