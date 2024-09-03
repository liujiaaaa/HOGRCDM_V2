FIT_HOGRCDM_Main<-function(ModelSetList,SizeList,Res,Q_H,Com_par=NULL){
  
 #########Check model setting 
  ModelSet<-
    ModelSetCheck(ModelSetList)
  
 # RegurlParaVec<-SettingList$RegurlParaVec
  
 ############The bottom layer
  
 if(ModelSetList$ModelDistri_Bottom=="Lognormal"){
   
   if( (ModelSetList$ModelFrame_Bottom=="Confirmatory"&ModelSetList$ModelFrame_Higher=="Confirmatory") ){
     
     FIT<-FIT_CC_Lognormal(Res, SizeList, InitAll, SettingList,
                                   ModelSetList, Q_H,Q_B)
     
   }else if((ModelSetList$ModelFrame_Bottom=="Exploratory" & ModelSetList$ModelFrame_Higher=="Exploratory")){
     
     
     
   }else if((ModelSetList$ModelFrame_Bottom=="Exploratory" & ModelSetList$ModelFrame_Higher=="Confirmatory")){
     
    
     
     
     RegurlParaVec<-SettingList$RegurlParaVec
     LL<-length(RegurlParaVec)
     skip_to_next<-rep(FALSE,LL)
     EBIC<-rep(NA,LL)
     
     
     for(lu in 1:LL){
       cat(paste("Runing algorithm with regularization parameter:",RegurlParaVec[lu],"\n","\n","\n"))
       skip_to_next[lu] <- FALSE
       
       # Note that print(b) fails since b doesn't exist
       RegurlPara<-rep(RegurlParaVec[lu],J);
       tryCatch({
         FIT<- FIT_CE_Lognormal(Res, SizeList, InitAll, SettingList,ModelSetList, Q_H,RegurlPara,Com_par)
         Par_size=sum(as.integer( FIT$Slope_B>RegurlParaVec[lu]/5))
         H_Num<-Par_size+sum(Q_Higher)+K;#The number of non-zero parameters
         EBIC[lu]<- -2*FIT$LIKELI+ H_Num*log(N)
         #+2*1*log(nchoosek_prac(sum(Q_Higher)+K+(K+1)*J,sum(Q_Higher)+K+Par_size))
       } , error = function(e) {
         cat("ERROR at Iteration", t,"","loop",lu)
         skip_to_next[lu]  <<- TRUE
       })
       
       if(skip_to_next[lu] ) { 
         EBIC[lu]<-Inf
         next 
       }
       
     }
     
     
     ll<-which.min(EBIC)
     RegurlPara<-rep(RegurlParaVec[ll],J);
     
     
     FIT<- FIT_CE_Lognormal(Res, SizeList, InitAll, SettingList,ModelSetList, Q_H,RegurlPara,Com_par)
     Par_size=sum(as.integer( FIT$Slope_B>RegurlParaVec[lu]/5))
     H_Num<-Par_size+sum(Q_Higher)+K;#The number of non-zero parameters
     EBIC_FIT<- -2*FIT$LIKELI+ H_Num*log(N) ###May be return this value?
     
     
   }
   
   
   
   
   
   
   
   
 }else if(ModelSetList$ModelDistri_Bottom=="LLM"){
   
   
   if( (ModelSetList$ModelFrame_Bottom=="Confirmatory"& ModelSetList$ModelFrame_Higher=="Confirmatory") ){
     
     FIT<-FIT_CC_LLM(Res, SizeList, InitAll, SettingList,
                           ModelSetList, Q_H,Q_B)
     
   }else if((ModelSetList$ModelFrame_Bottom=="Exploratory" & ModelSetList$ModelFrame_Higher=="Exploratory")){
     
     
     
   }else if((ModelSetList$ModelFrame_Bottom=="Exploratory" & ModelSetList$ModelFrame_Higher=="Confirmatory")){
     
     
     RegurlParaVec<-SettingList$RegurlParaVec
     LL<-length(RegurlParaVec)
     skip_to_next<-rep(FALSE,LL)
     EBIC<-rep(NA,LL)
     
     
     for(lu in 1:LL){
       cat(paste("Runing algorithm with regularization parameter:",RegurlParaVec[lu],"\n","\n","\n"))
       skip_to_next[lu] <- FALSE
       
       # Note that print(b) fails since b doesn't exist
       RegurlPara<-rep(RegurlParaVec[lu],J);
       tryCatch({
         FIT<- FIT_CE_LLM(Res, SizeList, InitAll, SettingList,ModelSetList, Q_H,RegurlPara,Com_par)
         Par_size=sum(as.integer( FIT$Slope_B>RegurlParaVec[lu]/5))
         H_Num<-Par_size+sum(Q_Higher)+K;#The number of non-zero parameters
         EBIC[lu]<- -2*FIT$LIKELI+ H_Num*log(N)
         #+2*1*log(nchoosek_prac(sum(Q_Higher)+K+(K+1)*J,sum(Q_Higher)+K+Par_size))
       } , error = function(e) {
         cat("ERROR at Iteration", t,"","loop",lu)
         skip_to_next[lu]  <<- TRUE
       })
       
       if(skip_to_next[lu] ) { 
         EBIC[lu]<-Inf
         next 
       }
       
     }
     
     
     ll<-which.min(EBIC)
     RegurlPara<-rep(RegurlParaVec[ll],J);
     
     
     FIT<- FIT_CE_LLM(Res, SizeList, InitAll, SettingList,ModelSetList, Q_H,RegurlPara,Com_par)
     Par_size=sum(as.integer( FIT$Slope_B>RegurlParaVec[lu]/5))
     H_Num<-Par_size+sum(Q_Higher)+K;#The number of non-zero parameters
     EBIC_FIT<- -2*FIT$LIKELI+ H_Num*log(N) ###May be return this value?
     
     
   }
     
     
   
   
   
   
   
   
   
 }else if(ModelSetList$ModelDistri_Bottom=="Poisson"){
   
   
   
   if( (ModelSetList$ModelFrame_Bottom=="Confirmatory"& ModelSetList$ModelFrame_Higher=="Confirmatory") ){
     
     FIT<-FIT_CC_Poisson(Res, SizeList, InitAll, SettingList,
                     ModelSetList, Q_H,Q_B)
     
   }else if((ModelSetList$ModelFrame_Bottom=="Exploratory" & ModelSetList$ModelFrame_Higher=="Exploratory")){
     
     
     
   }
 } else if(ModelSetList$ModelDistri_Bottom=="Gamma"){
     
     
     
     if( (ModelSetList$ModelFrame_Bottom=="Confirmatory"& ModelSetList$ModelFrame_Higher=="Confirmatory") ){
       
       FIT<-FIT_CC_Gamma(Res, SizeList, InitAll, SettingList,
                           ModelSetList, Q_H,Q_B)
       
     }else if((ModelSetList$ModelFrame_Bottom=="Exploratory" & ModelSetList$ModelFrame_Higher=="Exploratory")){
       
       
       
     }
   
   
   
   
   
   
   
   
   
   
   
   
   
 }
  
  
  
  
  
  
  return( FIT)
  
  
  
}
