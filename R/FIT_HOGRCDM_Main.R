FIT_HOGRCDM_Main<-function(ModelSetList,SizeList,Res,Q_H,Com_par=NULL,SettingList){
  
 # #########Check model setting 
 #  ModelSet<-
 #    ModelSetCheck(ModelSetList)
  
 # RegurlParaVec<-SettingList$RegurlParaVec
  
  if(ModelSetList$HigherLayer==T){
   FIT<- FIT_HOGRCDM_High_Bott(ModelSetList,SizeList,Res,Q_H,Com_par=NULL,SettingList)
  }else{
    FIT<- FIT_HOGRCDM_Bott(ModelSetList,SizeList,Res,Q_H,Com_par=NULL,SettingList)
  }
  return( FIT)
  
  
  
}
