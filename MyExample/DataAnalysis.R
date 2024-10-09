load("./LoadValue/DataAreBooklet1.RData")

Q_T<-Q_AreBooklet1
J<-nrow(Q_T)
Res1<-Res<-as.matrix(DataAreBooklet1[,-1])/60











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
  HigherLayer=F,
  ModelDistri_Bottom="Gamma",
  ModelStruc_Bottom="Main_effect",
  ModelStruc_Higher="Subscale",
  ModelFrame_Bottom="Confirmatory",
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
  HigherLayer=T,
  ModelDistri_Bottom=ModelSet$ModelDistri_Bottom,
  ModelStruc_Bottom=ModelSet$ModelStruc_Bottom,
  ModelStruc_Higher=ModelSet$ModelStruc_Higher,
  ModelFrame_Bottom=ModelSet$ModelFrame_Bottom,
  ModelFrame_Higher=ModelSet$ModelFrame_Higher
)

SizeList<-list(
  N=dim(Res1)[1],
  J=dim(Res1)[2],
  K=7,
  D=2
)



N<-SizeList$N
J<-SizeList$J
K<-SizeList$K
D<-SizeList$D
L<-2^K



Q_B<-Q_T
Q_H<-matrix(0,7,2)
Q_H[1:4,1]<-1
Q_H[5:7,2]<-1





TrueP_List_Bott<-NULL


####################This list is used for simulation with exploratory case where 
####################Hungarian algorithm is applied for computing Q_B matrix recovery
Com_par<-NULL








###Find initial values

InitAll<-INIT_HOGRCDM_Main(ModelSetList,SizeList,Res,Q_B,Q_H)
InitAll

InitAll$Slope_B_INI<-abs(InitAll$Slope_B_INI)
InitAll$Slope_H_INI<-abs(InitAll$Slope_H_INI)


SettingList<-list(
  StartSamSize=10,
  IncreSize=5,
  NumOfWeight=5,
  maxIter=100,
  ReturnLikelihood=T,
  epsCheck=0.04,
  Passing=3,
  PRINT=T,
  PR_ITER_NUM=5
)






###Fit Models

time_start<-proc.time()
FIT= FIT_HOGRCDM_Main(ModelSetList,SizeList,Res,Q_H,PenaltyList,SettingList)
time_end<-proc.time()-time_start


