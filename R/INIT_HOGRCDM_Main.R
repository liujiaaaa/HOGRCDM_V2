INIT_HOGRCDM_Main<-function(ModelSetList,SizeList,Res,Q_B,Q_H){
 if(ModelSetList$HigherLayer){
   InitialList<-INIT_HOGRCDM_Main_High_Bott(ModelSetList,SizeList,Res,Q_B,Q_H)
 }else{
   InitialList<- INIT_HOGRCDM_Main_Bott(ModelSetList,SizeList,Res,Q_B)
 }
  
  return(InitialList)
  
  
}