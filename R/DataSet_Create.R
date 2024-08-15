
DataSet_Create<-function(num_ob,num_attri,num_item,OriData){
  
  L<-2^num_attri
  Data_Layer<-matrix(NA,num_ob*L,num_item)
  for(l in 1:L){
    Data_Layer[ ((0:(num_ob-1))*L)+l,]<-OriData
  }
  return(Data_Layer)
}


# 
# Data_Layer0<-matrix(NA,N*L1,K0)
# for(l in 1:L1){
#   Data_Layer0[ ((0:(N-1))*L1)+l,]<-Res
# }
# 
# Data_Layer1<-matrix(NA,L1*L2,K1)
# for(l in 1:L2){
#   Data_Layer1[ ((0:(L1-1))*L2)+l,]<-Poss_Attr1
# }