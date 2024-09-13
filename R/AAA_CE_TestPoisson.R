


source_all <- function(directory = "R") {
  # Get a list of all R files in the specified directory
  files <- list.files(directory, pattern = "\\.R$", full.names = TRUE)
  
  # Source each file
  for (file in files) {
    source(file)
  }
}

# Example usage
source_all("R")



ModelSetList<-list(
  ModelDistri_Bottom="Poisson",
  ModelStruc_Bottom="Main_effect",
  ModelStruc_Higher="Subscale",
  ModelFrame_Bottom="Exploratory",
  ModelFrame_Higher="Confirmatory"
)

# ModelDistri_Bottom=ModelSet$ModelDistri_Bottom,
# ModelStruc_Bottom=ModelSet$ModelStruc_Bottom,
# ModelStruc_Higher=ModelSet$ModelStruc_Higher,
# ModelFrame_Bottom=ModelSet$ModelFrame_Bottom,
# ModelFrame_Higher=ModelSet$ModelFrame_Higher

LoadingPackages()

ModelSet<-
  ModelSetCheck(ModelSetList)

ModelSetList<-list(
  ModelDistri_Bottom=ModelSet$ModelDistri_Bottom,
  ModelStruc_Bottom=ModelSet$ModelStruc_Bottom,
  ModelStruc_Higher=ModelSet$ModelStruc_Higher,
  ModelFrame_Bottom=ModelSet$ModelFrame_Bottom,
  ModelFrame_Higher=ModelSet$ModelFrame_Higher
)

SizeList<-list(
  N=500,
  J=30,
  K=7,
  D=3
)



N<-SizeList$N
J<-SizeList$J
K<-SizeList$K
D<-SizeList$D
L<-2^K

######Loading Q matrices and True Parameters for Higher Order Structures: Slope_H and Interc_H
load("./LoadValue/QGU.RData")


if(ModelSetList$ModelStruc_Higher=="Subscale"){
  load("./LoadValue/Simple1.RData")
}else if(ModelSetList$ModelStruc_Higher=="Bifactor"){
  load("./LoadValue/BIFACTOR1.RData")
}


Q_B<-Q
Q_H<-t(Lam_1T)
Q_H[Q_H!=0]<-1


Slope_H<-Lam_1T
Interc_H<-Lam_0T

Sigma_thetaT<-diag(D)

Sigma_thetaT[2,3]<-Sigma_thetaT[3,2]<-0.35

TrueP_List_High<-list(
  Slope_H=Slope_H,
  Interc_H=Interc_H,
  Mu_thetaT=rep(0,D),
  Sigma_thetaT=Sigma_thetaT
  
)


##############This is only a procedure for my simulation
##############Users can input their own true parameters
##############This is only a procedure for my simulation
##############Users can input their own true parameters
if(ModelSetList$ModelStruc_Bottom=="Main_effect"){
  Beta1_MatT<-matrix(Q,J,K)
  Beta0T<-rep(NA,J)
  Beta_Para_List<-list()
  for(j in 1:J){
    Beta_Para_List[[j]]<-c(0.5, rep((1)/(rowSums(Q)[j]),rowSums(Q)[j]))
    Beta1_MatT[j,][Beta1_MatT[j,]==1]<- Beta_Para_List[[j]][-1]
    Beta0T[j]<-Beta_Para_List[[j]][1]
    #print(sum(Beta_Para_List[[j]]-c(Beta0T[j],Beta1_MatT[j,][Beta1_MatT[j,]!=0])))
  }
  Slope_B<-Beta1_MatT
  Interc_B<-Beta0T
  
  
}else if(ModelSetList$ModelStruc_Bottom=="All_effect"){
  
  
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
  
  
  # Beta0T<-rep(NA,J)
  SizeList$Q_Aug<-Q_Aug
  Beta1_MatT<-Q_Aug
  
  for( rr in 1:nrow( Beta1_MatT)){
    Beta1_MatT[rr,][ Beta1_MatT[rr,]!=0]<-1/(2^(sum(  Beta1_MatT[rr,(1:K)]))-1)
    # print(B_Slope1T[rr,][B_Slope1T[rr,]!=0])
  }
  Beta0T<-rep(0.5,J)
  
  Slope_B<-Beta1_MatT
  Interc_B<-Beta0T
  
  
  
  # 
  
  
}else if(ModelSetList$ModelStruc_Bottom=="DINA"){
  DINA_INDEX<-rep(1,J)
  Q_Aug<-Q
  q_num1<-K
  for(kk in 2:K){
    combinations <- t(combn(1:K, kk))
    te<-apply( Q,1,FUN= comFUN, combinations)
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
  
  
  # Beta0T<-rep(NA,J)
  SizeList$Q_Aug<-Q_Aug
  Beta1_MatT<-Q_Aug
  
  for( rr in 1:nrow( Beta1_MatT)){
    Beta1_MatT[rr,][ Beta1_MatT[rr,]!=0]<-1
    # print(B_Slope1T[rr,][B_Slope1T[rr,]!=0])
  }
  Beta0T<-rep(0.5,J)
  
  Slope_B<-Beta1_MatT
  Interc_B<-Beta0T
  
  
  
  
}



TrueP_List_Bott<-list(
  Slope_B=Slope_B,
  Interc_B=Interc_B,
  sd=rep(1,J),
  Q_Aug=Q_Aug
)


####################This list is used for simulation with exploratory case where 
####################Hungarian algorithm is applied for computing Q_B matrix recovery
Com_par<-list(
  Slope_B=TrueP_List_Bott$Slope_B,
  Interc_B=TrueP_List_Bott$Interc_B,
  SD=TrueP_List_Bott$sd,
  Slope_H=TrueP_List_High$Slope_H,
  Interc_H=TrueP_List_High$Interc_H,
  Mu_thetaT=TrueP_List_High$Mu_thetaT,
  Sigma_thetaT=TrueP_List_High$Sigma_thetaT
)


Res<-GenerateData(SizeList,ModelSetList,TrueP_List_Bott,TrueP_List_High,Q_B)








###Find initial values

InitAll<-INIT_HOGRCDM_Main(ModelSetList,SizeList,Res,Q_B=NULL,Q_H)
InitAll



SettingList<-list(
  StartSamSize=5,
  IncreSize=5,
  NumOfWeight=5,
  maxIter=100,
  ReturnLikelihood=T,
  epsCheck=0.04,
  Passing=3,
  PRINT=T,
  PR_ITER_NUM=5,
  RegurlParaVec=c(0.015,0.02)
)






###Fit Models
time_start<-proc.time()
FIT= FIT_HOGRCDM_Main(ModelSetList,SizeList,Res,Q_H,Com_par)
time_end<-proc.time()-time_start





#######Check the estimats
Test_Q<-FIT$Slope_B
Test_Q[Test_Q!=0]<-1
Q_Aug1<-TrueP_List_Bott$Slope_B
Q_Aug1[Q_Aug1!=0]<-1


####Should be 0
sum(Test_Q!=Q_Aug1)

####
ABias_SlopeB<-abs(FIT$Slope_B-TrueP_List_Bott$Slope_B)
ABias_InterB<-abs(FIT$Interc_B-TrueP_List_Bott$Interc_B)


ABias_SlopeH<-abs(FIT$Slope_H-TrueP_List_High$Slope_H)
ABias_InterH<-abs(FIT$Interc_H-TrueP_List_High$Interc_H)
#ABias_sd<-abs(FIT$SD-TrueP_List_Bott$sd)
ABias_Sigma<-abs(FIT$Sigma_theta-TrueP_List_High$Sigma_thetaT)



######Print results
mean(ABias_SlopeB[Q_Aug1!=0])
mean(ABias_SlopeB[Q_Aug1==0])#Should be 0
mean(ABias_InterB)
mean(ABias_SlopeH[t(Q_H)!=0])
mean(ABias_SlopeH[t(Q_H)==0])#Should be 0
mean(ABias_InterH)
#mean(ABias_sd)
sum(ABias_Sigma)/(SizeList$D^2-D)

print(time_end)

tt1<-which(colSums(TrueP_List_Bott$Slope_B)!=0)
tt2<-which(colSums(TrueP_List_Bott$Slope_B)==0)

library(pheatmap)
library(viridis)

par(mfrow = c(1, 2),  oma = c(0, 0, 0, 0))
heatmap(FIT$Slope_B[, tt1], Colv = NA, Rowv = NA, scale = "none", col = heat.colors(256))
heatmap(TrueP_List_Bott$Slope_B[, tt1], Colv = NA, Rowv = NA, scale = "none", col = heat.colors(256))


heatmap(FIT$Slope_B[, tt2], Colv = NA, Rowv = NA, scale = "none", col = heat.colors(256))
heatmap(TrueP_List_Bott$Slope_B[, tt2], Colv = NA, Rowv = NA, scale = "none", col = heat.colors(256))





