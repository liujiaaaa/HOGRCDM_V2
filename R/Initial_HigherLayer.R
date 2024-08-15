
Initial_HigherLayer<-function(Res, size_list,ModelSetList,
                              Beta0,Beta1_Mat0,Q_Ht){
  
  library(mvtnorm)
  
  model<-ModelSetList$ModelDistri_Bottom
  type<-ModelSetList$ModelStruc_Bottom
  
  
  if(ModelSetList$ModelStruc_Higher=="Bifactor"){
    isBifactor<-T
  }else{
    isBifactor<-F
  }
  
  D=size_list$D
  
  
  # if(isBifactor){
  #   D=size_list$D-1
  # }else{
  #     D=size_list$D
  #   }
  
  time_start<-proc.time()
  library(mvnfast)
  
  N<-size_list$N
  J<-size_list$J
  K<-size_list$K
  # D<-size_list$D
  #maxIter<-size_list$maxIter
  L<-2^K
  
  
  # Alpha<-Initial_List$Alpha0
  #theta<-Initial_List$theta0[1:L,]
  # Lam_1<-Initial_List$Lam_1_0
  # Lam_0<-Initial_List$Lam_0_0
  
  
  
  
  
  Beta0<-Beta0
  Beta1_Mat<-Beta1_Mat0
  
  
  if(type=="GDINA"|type=="DINA"){
    Q_Aug<-size_list$Q_Aug
    InitListUpdate<-list()
    InitListUpdate$Beta1_Mat<-matrix(0,nrow(Q_Aug),ncol(Q_Aug))
    InitListUpdate$Beta1_Mat[,1:K]<-Beta1_Mat
    InitListUpdate$Beta0<-Beta0
    
    
    Beta0<-InitListUpdate$Beta0
    Beta1_Mat<-InitListUpdate$Beta1_Mat
  }
  
  
  
  
  #Beta1_Mat<-Beta1_Mat*Q
  # 
  # Lam_1<-Lam_1T
  # Lam_0<-Lam_0T
  
  # Beta0<-Beta0T
  # Beta1_Mat<-Beta1_MatT
  # 
  #PPL<- Initial_List$PPL
  
  #isBifactor<-SettingList$isBifactor
  # Pre_period<-SettingList$Pre_period
  # Subsize<-SettingList$Subsize
  # maxIter<-SettingList$maxIter
  # burnin<-SettingList$burnin
  # gamma<-SettingList$gamma
  # LLambda<-SettingList$LLambda
  # Likelihood<-SettingList$Likelihood
  # 
  
  
  
  index<-list()
  for(k in 1:K){
    index[[k]]<- which(Q_Higher[k,]==1)
  }
  
  #L<-2^K
  Poss_Attr<-matrix(NA,L,K)
  for(l in 1:K){
    repTimes<-2^(l-1)
    Poss_Attr[,l]<-rep(c(rep(1,L/2^l),rep(0,L/2^l)),repTimes)
  }
  
  NPoss_Attr<-1-Poss_Attr
  
  
  
  
  if(type=="Main_effect"){
    Data_x1<- matrix(NA,K,L*N)
    Data_x1[,]<-t(Poss_Attr)
    Data_x<-t(Data_x1)
    
    Data_Res<-matrix(NA,N*L,J)
    for(l in 1:L){
      Data_Res[ ((0:(N-1))*L)+l,]<-Res
    }
    Poss_AttrUp<- Poss_Attr
    
  }else if(type=="GDINA"|type=="DINA"){
    Poss_Attr1<-PossATTRIB(K)
    ####################################################
    PS<-AugmenInteraction(OrigriData=Poss_Attr1,K=K)
    Poss_Attr1_Aug<-PS$OrigriData_Aug
    s_num1<-PS$s_num
    
    
    Data_x<-Data_Pred(num_ob=N,num_attri= K,num_attri_aug= s_num1,
                      Ori_Attri= Poss_Attr1_Aug)
    ############################################################
    
    Data_Res<-DataSet_Create(num_ob=N,num_attri=K,num_item=J,OriData=Res)
    
    Poss_AttrUp<- Poss_Attr1_Aug
  }
  
  
  
  
  
  
  
  
  
  
  IA2<- array((-1)^t((Poss_Attr[, rep(1:K, each = K)] +
                        Poss_Attr[, rep(1:K, times = K)])), dim = c(K, K, L))
  
  
  IA<-(-1)^NPoss_Attr
  
  
  #ASA<-array(NA,dim=c(K,K,L))
  #TsASA<-array(NA,dim=c(K,K,L))
  
  ST<-array(NA,dim=c(K,K,L))
  
  IAd<-matrix(NA,L,K)
  SdASM<-matrix(NA,L,K)
  PA<-rep(NA,L)
  LOW<-rep(-Inf,K)
  MEAN<-rep(0,K)
  Beta0_matrix<-matrix(NA,L,J)
  PRA<-array(NA,dim=c(L,J,N))
  PR1A<-PR2A<-array(NA,dim=c(L,J,N))
  ResA<-array(NA,dim=c(L,J,N))
  Lam_0_matrix<-matrix(NA,L,K)
  for(l in 1:L){
    ResA[l,,]<-t(Res)
  }
  ResAN<-1-ResA
  
  ResM<-Res[rep(1:nrow(Res),L),]
  rLN<-rep(1:L,each=N)
  ResMN<-1-ResM
  
  PAM<-matrix(NA,L,N)
  
  
  
  
  
  I1L<-rep(1,L)
  I1N<-rep(1,N)
  PRTNJL<-matrix(NA,N,L)
  
  
  
  PRT1<-matrix(0,L*N,J)
  
  # num_cores <- detectCores()-1
  #  cl <- makeCluster(num_cores)
  # registerDoParallel(cl)
  
  
  
  
  #if(i%%10==0){
  
  # print(Beta0)
  # print( Beta1_Mat)
  # print( Lam_1)
  # print( Lam_0)
  # 
  # #print(ta)
  # print(Sigma_theta)
  #  print(Lambda)
  # }
  
  
  
  
  
  Beta0_matrix[,]<-I1L%*%t(Beta0)
  t_Beta1_Mat<-t(Beta1_Mat)
  
  AB<-Poss_AttrUp%*%t_Beta1_Mat+Beta0_matrix
  
  if(model=="Poission"|model=="Gamma"){
    #########################
    AB[AB<0]<- 0.1
    ##########################
  }
  
  if(model=="Poission"){
    PRT1[,]<-dpois(as.vector(ResM),as.vector(AB[rLN,]))
  }else if(model=="LLM"){
    
    Beta0_matrix[,]<-I1L%*%t(Beta0)
    t_Beta1_Mat<-t(Beta1_Mat)
    PR1<-plogis(AB)
    
    PRt1<- PR1[rLN,]
    PRt2<-1-PRt1
    PRT1<-ResM*PRt1+ResMN*PRt2
  }else if(model=="Lognormal"){
    sdd<-apply(log(ResM),2,sd)
    PRT1[,]<-dnorm(as.vector(log(ResM)),as.vector(AB[rLN,]),sd= sdd)
  }else if(model=="Gamma"){
    
    mm<-apply(Res,2,mean)
    
    ss<-apply(Res,2,var)
    
    iniShape<-mm^2/ss
    
    
    
    PRT1[,]<-dgamma(as.vector(ResM),shape= iniShape,rate=AB[rLN,])
    
  }
  
  
  
  
  
  
  
  
  
  
  # PR2<-1-PR1
  #  PR1A[,,]<-PR1
  #  PR2A[,,]<-PR2
  # PRT<-ResA*PR1A+ResAN*PR2A
  #PRTNJ<-apply(PRT,c(1,3),prod)
  
  
  # PRt1<- PR1[rLN,]
  # PRt2<-1-PRt1
  #PRT1<-ResM*PRt1+ResMN*PRt2
  #PRTNJ<-matrix(rowProds(PRT1),L,N,byrow=T)
  PRTNJL[,]<-rowProds(PRT1)
  PRTNJ<-t(PRTNJL)
  
  
  
  
  NUM<-PRTNJ
  
  DEM<-colSums(NUM)
  DEMM<-rep(1,L)%*%t(DEM)
  
  
  PRAA<-NUM/DEMM 
  #default family is gaussian
  ppA<-(rowMeans(PRAA))
  
  
  samTh<-sample(1:L,size=30000,replace=T,prob=ppA)
  
  DATA_ALPHA<-Poss_Attr[samTh,]
  
  
  InitL<- IniFUN_Probit(DATA=DATA_ALPHA, N=length(samTh),J=K,K=D,eplison_Layer0=1e-04,TARGET=TARGET)
  
  Lam_1_0<-t(InitL$Init_B_Slope1)
  Lam_0_0<-InitL$Init_B_Intec1
  
  
  
  
  Lam_1_te<-abs(Lam_1_0)
  for(t in 1: ncol(Lam_1_te)){
    Lam_1_te[,t]<-Lam_1_te[,t]/max(Lam_1_te[,t])
  }
  Lam_1_teO<-Lam_1_te
  
  compute_cost_matrix <- function(Q_Ht, Lam_1_te) {
    n <- nrow(Q_Ht)
    cost_matrix <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        cost_matrix[i, j] <- sum((Q_Ht[i, ] - Lam_1_te[j, ])^2)
      }
    }
    cost_matrix
  }
  
  
  # Initial column order
  initial_col_order <- 1:ncol(Lam_1_te)
  initial_row_order <- 1:nrow(Lam_1_te)
  
  # Iterative refinement
  max_iter <- 100
  tol <- 1e-6
  for (iter in 1:max_iter) {
    # Compute row permutation
    row_cost_matrix <- compute_cost_matrix(Q_Ht, Lam_1_te)
    row_perm <- solve_LSAP(row_cost_matrix)
    Lam_1_te <- Lam_1_te[row_perm, ]
    
    # Compute column permutation
    col_cost_matrix <- compute_cost_matrix(t(Q_Ht), t(Lam_1_te))
    col_perm <- solve_LSAP(col_cost_matrix)
    Lam_1_te <- Lam_1_te[, col_perm]
    
    # Update column order
    initial_col_order <- initial_col_order[col_perm]
    initial_row_order <- initial_row_order[row_perm]
    
    # Check for convergence (optional)
    if (max(abs(row_perm - 1:nrow(Q_Ht)), abs(col_perm - 1:ncol(Q_Ht))) < tol) {
      break
    }
  }
  
  
  
  
  Lam_1_0<- Lam_1_0[, initial_col_order][ initial_row_order,]
  Lam_0_0<- Lam_0_0[initial_col_order]
  # Print the reordered matrix Lam_1_te
  #print(Lam_1_te)
  
  
  
  
  
  return(list(Lam_1= Lam_1_0,Lam_0= Lam_0_0,  DATA_ALPHA=  DATA_ALPHA))
  
  
}
