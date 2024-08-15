Data_Pred<-function(num_ob,num_attri,num_attri_aug,Ori_Attri){
  L<-2^num_attri
  Data_x1<- matrix(NA,num_attri_aug,L*num_ob)
  Data_x1[,]<-t(Ori_Attri)
  Data_x<-t(Data_x1)
  return(Data_x)
}
