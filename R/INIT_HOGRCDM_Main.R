INIT_HOGRCDM_Main<-function(ModelSetList,SizeList,Res,Q_B,Q_H){
  N<-SizeList$N
  J<-SizeList$J
  K<-SizeList$K
  D<-SizeList$D
  L<-2^K
  
  Q_Ht<-t(Q_H)
  
  if(ModelSetList$ModelStruc_Higher=="Bifactor"){
    isBifactor<-T
  }else{
    isBifactor<-F
  }
  
  
if(ModelSetList$ModelDistri_Bottom=="Lognormal"){
  

    if( (ModelSetList$ModelFrame_Bottom=="Confirmatory" & ModelSetList$ModelFrame_Higher=="Confirmatory") ){
      
      
      ########################Initilization for bottom layer
      INI_B<-INI_CC_Lognormal(DATA=log(Res), N,J,K,eplison_Layer0=1e-04)
      
      Slope_B_INI<-INI_B$Init_B_Slope1
      Interc_B_INI<-INI_B$Init_B_Intec1
      # Interc_B_INI[Interc_B_INI<0]<- 0
      # Slope_B_INI[ Slope_B_INI<0]<-0
      Q_B_INI<- Slope_B_INI
      Q_B_INI[ Q_B_INI!=0]<-1
      
      
      
      KK=ncol( Q_B)
      DD<-nrow( Q_B)
     
      Cost_matrix<-matrix(NA,KK,KK)
      for(ll in 1:KK){
        for(mm in 1:KK){
          Cost_matrix[ll,mm]<-sum(abs(Q_B_INI[,mm]-Q_B[,ll]))
        } 
      }
      Test<-HungarianSolver(Cost_matrix) 
      if(sum(abs(Test$pairs[,2]-(1:KK)))!=0) {
        print(Test$pairs[,2]);Hungarian=T
      }else{
        Hungarian=F
      }
      Slope_B_INI<-  Slope_B_INI[, Test$pairs[,2]]
      
      
      ########################Initilization for higher order layer   
      
      INI_H<-INI_CC_HigherLayer(Res, size_list=SizeList,ModelSetList,
                                     Beta0=Interc_B_INI, Beta1_Mat0=Slope_B_INI,Q_Ht=t(Q_H),Ini_AlphaSize=5000)
      
      
      Slope_H_INI=INI_H$Lam_1 ####May be a target rotation?
      Interc_H_INI=INI_H$Lam_0
      
      
      
      
      SD<-apply(log(Res),2,sd)
      
      
     
      
      
      
      
      
      
    }else{
      
      
      
      
    }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
  

  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  }else if(ModelSetList$ModelDistri_Bottom=="LLM"){
 
  
  
  
  SD=NULL 
}else if(ModelSetList$ModelDistri_Bottom=="Poisson"){
  
  
  
  
  SD=NULL
}
  
  
  
  
  
  
  
  
  
  
  if((ModelSetList$ModelStruc_Bottom=="All_effect")){
    
    Slope_B_INI_Aug<- matrix(0, J, L-1)
    Slope_B_INI_Aug[,1:K]<- Slope_B_INI
    Slope_B_INI<-  Slope_B_INI_Aug
    
  }else if(ModelSetList$ModelStruc_Bottom=="DINA"){
    
    ##################May be  I should adjust the initilization scheme for DINA case?
    
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
    
    
    
    
    Slope_B_INI_Aug<- matrix(0, J, L-1)
   for(jj in 1:J){
     Slope_B_INI_Aug[jj,Q_Aug[jj,]!=0]<-sum( Slope_B_INI[jj,])
   }
    Slope_B_INI<-  Slope_B_INI_Aug
  }
  
  
  
  
  
  InitialList<-list(
    Slope_H_INI=Slope_H_INI,
    Interc_H_INI=Interc_H_INI,
    Slope_B_INI=Slope_B_INI,
    Interc_B_INI=Interc_B_INI,
    Sigma_theta_INI=diag(D),
    SD=SD
  )
  
  return(InitialList)
  
  
}