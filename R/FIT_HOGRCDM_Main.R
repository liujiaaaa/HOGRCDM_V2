FIT_HOGRCDM_Main<-function(ModelSetList,SizeList,Res,Q_H,PenaltyList){
  
 #########Check model setting 
  ModelSet<-
    ModelSetCheck(ModelDistri_Bottom=ModelSetList$ModelDistri_Bottom,
                  ModelStruc_Bottom=ModelSetList$ModelStruc_Bottom, ModelStruc_Higher=ModelSetList$ModelStruc_Higher)
  
  
 ############The bottom layer
  
 if(ModelSetList$ModelDistri_Bottom=="Lognormal"){
   
   if( (ModelSetList$ModelFrame_Bottom=="Confirmatory"&ModelSetList$ModelFrame_Higher=="Confirmatory") ){
     
     FIT<-FIT_CC_Lognormal(Res, SizeList, InitAll, SettingList,
                                   ModelSetList, Q_H,Q_B)
     
   }else{
     
   }
   
   
   
   



   
   
   
   
   
 }else if(ModelSetList$ModelDistri_Bottom=="LLM"){
   
   
   
   
 }else if(ModelSetList$ModelDistri_Bottom=="Poisson"){
   
   
   
   
 }
  
  
  
  
  
  
  
  
  
}
