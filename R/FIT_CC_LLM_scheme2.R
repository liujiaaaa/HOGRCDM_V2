
FIT_CC_LLM_scheme2<-function(Res, SizeList, InitAll, SettingList,
                           ModelSetList, Q_H, RegurlPara
){
  
  time_start<-proc.time() 
  
  
  # Opt_Weight_Higher<-function(j,x,Q_mat,Res,LLambda,Link="identity", weights){
  #   data_x<-x[,which(Q_mat[j,]!=0)]
  #   fit <- glmnet(x=data_x,Res[,j],  family = binomial(link = Link), weights=weights,
  #                 lambda=LLambda[j])
  #   as.vector(coef(fit))
  # }
  
  
  
  
  # E-step to obtain mc_samples of factors for each observation
  sample_single_thetas <- function(m, k, mc_samples, d1, d2, s_arr,Sigma_theta){
    custom_library_path <- "/burg/home/jl6795/R/x86_64-pc-linux-gnu-library/4.3"
    library(TruncatedNormal,lib.loc = custom_library_path)
    inv_term <- solve(d1%*%Sigma_theta%*%t(d1)+diag(m))
    cov_vo <- Sigma_theta - Sigma_theta%*%t(d1)%*%inv_term %*% d1%*%Sigma_theta
    cov_v1 <- diag(1/s_arr) %*% (d1%*%Sigma_theta%*%t(d1) + diag(m)) %*% diag(1/s_arr)
    trunc_lev <- -diag(1/s_arr) %*% d2
    v0_s <- MASS::mvrnorm(n=mc_samples, mu = rep(0, k), cov_vo)
    v1_s <-TruncatedNormal::rtmvnorm(n=mc_samples, mu = rep(0, m), sigma= cov_v1, lb = trunc_lev )
    linear_term <- t(d1)%*%inv_term%*%diag(s_arr)
    theta_i <- v0_s +t(linear_term %*% t(v1_s))%*%Sigma_theta
    return(theta_i)
  }
  
  
  
  
  
  
  
  sample_thetas <- function(L,K,D,data,Slope_H,Interc_H,Sigma_theta,  mc_samples){
    #source("sample_single_thetas.R")
    # source("/burg/home/jl6795/DINA2/DINA/sample_single_thetas.R")
    n = L
    m =K
    k =D
    library(magrittr)
    
    d1_array <- sweep((rep(1, n) %x% t( Slope_H)), 1, c(t(2*data-1)), "*")
    d2_array <- sweep(2*data-1, 2, (Interc_H), "*")
    
    s_array<- sqrt(diag(t(Slope_H)%*%Sigma_theta%*%Slope_H)+1)
    
    # s_array <- (rowSums(t(Slope_H)^2)+1)^0.5
    thetas <- 1:n %>% 
      furrr::future_map(function(.x) sample_single_thetas(m, k, mc_samples, d1_array[((.x-1)*m+1):(.x*m), ], d2_array[.x, ], s_array,Sigma_theta), .options = furrr::furrr_options(seed = TRUE)) %>% 
      do.call(rbind, .)
    return(thetas)
  }
  
  library(data.table)
  
  FindList <- function(j, data_y, data_x, Q_Bm) {
    index_j <- which(Q_Bm[j,] != 0)
    A <- as.matrix(data_y[,j])  # Response
    B <- as.matrix(data_x[, index_j])  # Attributes
    AB<-cbind(A,B)
    
    # Combine A and B into one data.table
    dt_AB <- data.table(AB)
    
    # Create a combined key from both A and B for identifying identical rows
    dt_AB[, key := apply(.SD, 1, paste, collapse = "_"), .SDcols = names(dt_AB)]
    
    data_u1<-unique(dt_AB$key)
    
    # Get unique rows and their corresponding indices
   dt_unique <- unique(dt_AB, by = "key")
   dt_uniqueM<-as.matrix(dt_unique)
   dt_uniqueM<-dt_uniqueM[,-ncol(dt_uniqueM)]
   dt_uniqueM<- matrix(as.numeric( dt_uniqueM),ncol=ncol(AB),byrow=F)
    # Map the combined_key back to its unique index group
    INDEX_Wei <- match(dt_AB$key, dt_unique$key)
    
    
    unique_values <- unique(INDEX_Wei)
    # Find the indices for each unique value
    
    
    INDEX_Wei1<-INDEX_Wei
    for(xx in  unique_values ){
      INDEX_Wei1[ which(INDEX_Wei == xx)[-1]]<-0
    }
    

    UpA<- dt_uniqueM[,1:(ncol(AB)-1)]
    UpB<- dt_uniqueM[,ncol(AB)]
   
    
 #indices<-which(INDEX_Wei1!=0)
    
    return(list(INDEX_Wei = INDEX_Wei, dRes = UpA, dX = UpB))
  }
  

  
  ######################## 
  #lRes<-log(Res)
  ######################## 
  
  N<-SizeList$N
  J<-SizeList$J
  K<-SizeList$K
  D<-SizeList$D
  L<-2^K
  
  NumOfWeight=SettingList$NumOfWeight
  StartSamSize=SettingList$StartSamSize
  IncreSize=SettingList$IncreSize
  NumOfWeight=SettingList$NumOfWeight
  ModelStruc_Bottom=SettingList$ModelStruc_Bottom
  isBifactor=ModelSetList$ModelStruc_Higher=="Bifactor"
  maxIter=SettingList$maxIter
  ReturnLikelihood=SettingList$ReturnLikelihood
  epsCheck=SettingList$epsCheck
  Passing=SettingList$Passing
  ModelStruc_Bottom=ModelSetList$ModelStruc_Bottom
  ModelDistri_Bottom=ModelSetList$ModelDistri_Bottom
  ModelStruc_Higher=ModelSetList$ModelStruc_Higher
  
  
  
  
  #########Initilization
  Slope_H<-InitAll$Slope_H_INI
  Interc_H<-InitAll$Interc_H_INI  
  Interc_B<-InitAll$Interc_B_INI 
  Slope_B<-InitAll$Slope_B_INI
  Sigma_theta<-InitAll$Sigma_theta_INI
  ta<-InitAll$SD#initial values for sd
  
  
  
  if(SettingList$PRINT==T){
    PRINT=T
    
    if(is.null(SettingList$PR_ITER_NUM)){
      PR_ITER_NUM=5
    } else{
      PR_ITER_NUM=SettingList$PR_ITER_NUM
    }
    
  } else{
    PRINT=F
  }
  
  
  
  ###########################Testing
  # Slope_B<-Com_par$Slope_B
  # Interc_B<-Com_par$Interc_B
  # Slope_H<-Com_par$Slope_H
  # Interc_H<-Com_par$Interc_H
  
  #L<-2^K
  Poss_Attr<-PossATTRIB(K)
  
  NPoss_Attr<-1-Poss_Attr
  
  
  
  
  if(ModelStruc_Bottom=="Main_effect"){
    Q_Aug<-Q_B
    
    Poss_AttrUp<- Poss_Attr
    
    Data_x1<- matrix(NA,K,L*N)
    Data_x1[,]<-t(Poss_AttrUp)
    Data_x<-t(Data_x1)
    
    Data_Res<-matrix(NA,N*L,J)
    for(l in 1:L){
      Data_Res[ ((0:(N-1))*L)+l,]<-Res
    }
    
    
    
  }else  if(ModelStruc_Bottom=="All_effect"){
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
    
    
  }else if (ModelStruc_Bottom=="DINA"){
    
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
  
  
  
  
  
  
  
  
  
  DATA_INDEX_LIST<-lapply(1:J,FindList,data_y=Data_Res, data_x=Data_x,Q_Bm=Q_B)
  
  
  
  
  #Slope_H<-Slope_HT
  #Interc_H<-Interc_HT
  # Slope_B<-Beta1_MatT
  #Interc_B<-Beta0T
  
  index<-list()
  for(k in 1:K){
    index[[k]]<- which(Q_H[k,]==1)
  }
  
  
  #Data_x<-Poss_Attr
  #for(su in 2:N){
  #  Data_x<-rbind(Data_x,Poss_Attr)
  #}
  
  
  
  
  LAMBDA<-rbind(Interc_H,Slope_H)
  LAMBDA[-1,][t(Q_H)==0]<-0
  QQc<-cbind(Q_H,rep(1,K))
  
  
  
  
  T7<-0
  T1<-T2<-0
  
  IA2<- array((-1)^t((NPoss_Attr[, rep(1:K, each = K)] +
                        NPoss_Attr[, rep(1:K, times = K)])), dim = c(K, K, L))
  
  
  IA<-(-1)^NPoss_Attr
  
  
  ST<-array(NA,dim=c(K,K,L))
  
  IAd<-matrix(NA,L,K)
  SdASM<-matrix(NA,L,K)
  PA<-rep(NA,L)
  LOW<-rep(-Inf,K)
  MEAN<-rep(0,K)
  Interc_B_matrix<-matrix(NA,L,J)
  PRA<-array(NA,dim=c(L,J,N))
  PR1A<-PR2A<-array(NA,dim=c(L,J,N))
  ResA<-array(NA,dim=c(L,J,N))
  
  Interc_H_matrix<-matrix(NA,L,K)
  # for(l in 1:L){
  #   ResA[l,,]<-t(Res)
  # }
  # ResAN<-1-ResA
  
  ResM<-Res[rep(1:nrow(Res),L),]
  ResMN<-1-ResM
  rLN<-rep(1:L,each=N)
  #repmat(t(Res), 1,L)
  #replicate(L,t(Res))
  #  ResMN<-1-ResM
  
  PAM<-matrix(NA,L,N)
  
  
  
  
  I1L<-rep(1,L)
  I1N<-rep(1,N)
  
  PRTNJL<-matrix(NA,N,L)
  
  
  
  PRT1<-matrix(0,L*N,J)
  
  PACDF<-function(ll,LOW,upperMat,corrArr){
    mvtnorm::pmvnorm(lower=LOW,upper=upperMat[ll,],mean=MEAN,corr=corrArr[,,ll])
  }
  
  # num_cores <- detectCores()-1
  #  cl <- makeCluster(num_cores)
  # registerDoParallel(cl)
  
  eps=1
  i=1
  BetaP<-NULL
  
  ##############The first Pre_period iterations,  do not specify warm start for 
  ##############optimization
  
  # cl <- makeCluster(detectCores() - 1)  # Reserve one core for system processes
  # registerDoParallel(cl)
  # 
  
  Beta<-rbind(Interc_B,t(Slope_B) )
  
  passing=0
  
  while((passing<Passing & i<maxIter)){
    
    if((PRINT&i%%PR_ITER_NUM==0)){
      print(Slope_B)###This can be revised to print more information
      print(Interc_B)
      print( Slope_H)
      print( Interc_H)
      print(Sigma_theta)
      print(ta)
    }
    
    # print(Slope_B)###This can be revised to print more information
    # print(Interc_B)
    # print( Slope_H)
    # print( Interc_H)
    # print(Sigma_theta)
    # print(ta)

    
   print(eps)
    
    
    Pre_LAMBDA<-LAMBDA
    Pre_Sigma_theta<-Sigma_theta
    Pre_Beta<-  Beta
    
    AS<-t(Slope_H)%*%Sigma_theta%*%Slope_H
    
    SdAS<-sqrt(1+diag(AS))
    SdASM[,]<-I1L%*%t(SdAS)
    TsAS<-SdAS%*%t(SdAS)
    # TsASA[,,]<-TsAS
    
    AT<-AS/TsAS
    diag( AT)<-1 #cor of diagonal equal to 1
    
    ST[,,]<- AT
    #SA is covariance matrix
    SA<-IA2*ST
    
    #TsASA[,,]<-TsAS
    # ASA[,,]<-AS
    # SA<-(ASA*IA2)/TsASA
    
    #for(l in 1:L){
    #  diag(SA[,,l])<-1
    #}
    
    IAd<-IA*(I1L%*%t(Interc_H))
    Fq<-IAd/SdASM
    
    # for(l in 1:L){
    #   PA[l]<- mvtnorm::pmvnorm(lower=LOW,upper=Fq[l,],mean=MEAN,corr=SA[,,l])
    # }
    
    
    
    PA<-apply(matrix(1:L),1,FUN=PACDF,LOW=LOW,upperMat=Fq,corrArr=SA)  
    
    ###################################################  
    if(i<NumOfWeight)  PA<-PA*1/2+rep(1/L,L)*1/2
    ###################################################  
    
    
    
    #  PA <- foreach(l = 1:L, .combine = "c") %dopar% {
    #   library(mvtnorm)
    #   pmvnorm(lower = LOW, upper = Fq[l,], mean = MEAN, corr = SA[,,l])
    # }
    
    
    PAM[,]<-PA%*%t(I1N)
    
    
    Interc_B_matrix[,]<-I1L%*%t(Interc_B)
    t_Beta1_Mat<-t(Slope_B)
    
    AB<-Poss_AttrUp%*%t_Beta1_Mat+Interc_B_matrix
    #########################
    #  AB[AB<0]<- 0
    ##########################
    
    PR1<-plogis(Poss_AttrUp%*%t_Beta1_Mat+Interc_B_matrix)
    
    PRt1<- PR1[rLN,]
    PRt2<-1-PRt1
    PRT1<-ResM*PRt1+ResMN*PRt2
  
    PRTNJL[,]<-rowProds(PRT1)
    PRTNJ<-t(PRTNJL)
    
    
    
    
    NUM<-PRTNJ*PAM
    
    DEM<-colSums(NUM)
    DEMM<-rep(1,L)%*%t(DEM)
    
    #Posterior
    PRAA<-NUM/DEMM 
    #default family is gaussian
    
    
    
    
    
    ppA<-(rowSums(PRAA))
    
    PRAV<-sqrt(ppA)
    PRAVM<-PRAV%*%t(rep(1,D))
    WEI<-as.vector(PRAA)
    
    PRAV_s<-as.vector(t(PRAA))%*%t(rep(1,J))
    
    
    
# 
    Beta<-sapply(1:J, OptBott_ConfirUp,  DATA_INDEX_LIST,Q_mat=Q_Aug,weights=WEI, N,Family="binomial", Link="logit",
                 simplify = T)


   
    # 
    # Beta<-sapply(1:J, OptBott_Confir, x=Data_x,
    #              Res=Data_Res,Q_mat=Q_Aug,weights=WEI,Family="binomial", Link="logit",
    #              simplify = T)

    Interc_B <- Beta[1,]
    Slope_B <- t(Beta[-1,])
    
    
    
    
    
    
    # 
    # time1 <- system.time({
    #   
    # })
    
    
    # time2 <- system.time({
    #   Beta1<-sapply(1:J, OptBott_ConfirWithMerge, x=Data_x,
    #                 ResD=Data_Res,Q_mat=Q_Aug,weights=WEI,Family="binomial", Link="logit",
    #                 simplify = T)
    #   Interc_B1 <- Beta1[1,]
    #   Slope_B1 <- t(Beta1[-1,])
    # })
   
    
    
 
   
    
    
    
    colnames(Interc_B)<-NULL
    rownames(Slope_B)<-colnames( Slope_B)<-NULL
    
    
    
    
    thetas1<- sample_thetas(L,K,D,data=Poss_Attr,Slope_H,Interc_H,
                            Sigma_theta=Sigma_theta,mc_samples=IncreSize)
    thetaA<-array(NA,dim=c(L,IncreSize,D))
    for(ll in 1:L){
      # print((((ll-1)*mc_samples)+1):(ll*mc_samples))
      thetaA[ll,,]<-thetas1[(((ll-1)*IncreSize)+1):(ll*IncreSize),]
    }
    thetaA1<-aperm( thetaA,c(3,2,1))
    
    
    theta_At<-matrix(NA,D,L*IncreSize)
    theta_At[,]<-thetaA1
    thetaD<-t(theta_At[,])
    
    xxD<-Poss_Attr[rep(1:L,each=IncreSize),]
    
    wwe<-rep(ppA,each=IncreSize)   
    
    wee1<- rep(PRAV,each=IncreSize)  
    
    thetas2<- thetas1*(wee1%*%t(rep(1,D)))
    
    
    
    
    
    # LAMBDA<-sapply(1:K,  Opt_Weight_Higher,Q_mat=Q_H,
    #                x=thetaD,
    #                Res=xxD,weights= wwe, Link="probit",
    #                LLambda=rep(0,K), simplify = T)
    
    LAMBDA<-sapply(1:K,    Opt_Weight_Higher_Simpler,Q_mat=Q_H,
                   x=thetaD,
                   Res=xxD,weights= wwe, Link="probit",
                   simplify = T)
    
    
    
    #print(test)
    
    # result_Lamnda<-rbind(LAMBDA[-1,],LAMBDA[1,])
    #  Lambda[t(QQc)!=0]<-result_Lamnda
    #  print( Lambda)
    
    
    Slope_H<-LAMBDA[-1,]
    Interc_H<-LAMBDA[1,]  
    
    # Sigma_theta<-(theta_At)%*%t(theta_At)/(L*IncreSize)
    Sigma_theta<-t(thetas2)%*%thetas2/(IncreSize*N)
    
    if(isBifactor){
      Sigma_theta1<-Sigma_theta[-1,-1]
      D1<-D-1
      Zrs<-tanh(Sigma_theta1)
      diag(Zrs)<-0
      Zrs[lower.tri(Zrs)]<-0
      
      UU<-matrix(NA,D1,D1)
      UU[lower.tri(UU)]<-0
      UU[1,1]<-1
      UU[1,-1]<-Zrs[1,-1]
      
      for(r in 2:D1){
        for(s in r:D1){
          
          if(r==s){
            UU[r,s]<-prod((1-Zrs[(1:(r-1)),s]^2)^(1/2))
          } else{
            UU[r,s]<- Zrs[r,s]* prod((1-Zrs[(1:(r-1)),s]^2)^(1/2))
          }
          
          
        }
      }
      Sigma_theta[-1,-1]<-t(UU)%*%UU
      Sigma_theta[1,-1]<-0
      Sigma_theta[-1,1]<-0
      Sigma_theta[1,1]<-1
      
      # Sigma_theta[-1,-1]<-cov2cor(Sigma_theta[-1,-1])
      # Sigma_theta[1,-1]<-0
      # Sigma_theta[-1,1]<-0
      # Sigma_theta[1,1]<-1
    }else
    {
      
      #################STAN 1
      Sigma_theta1<-Sigma_theta
      D1<-D
      Zrs<-tanh(Sigma_theta1)
      diag(Zrs)<-0
      Zrs[lower.tri(Zrs)]<-0
      
      UU<-matrix(NA,D1,D1)
      UU[lower.tri(UU)]<-0
      UU[1,1]<-1
      UU[1,-1]<-Zrs[1,-1]
      
      for(r in 2:D1){
        for(s in r:D1){
          
          if(r==s){
            UU[r,s]<-prod((1-Zrs[(1:(r-1)),s]^2)^(1/2))
          } else{
            UU[r,s]<- Zrs[r,s]* prod((1-Zrs[(1:(r-1)),s]^2)^(1/2))
          }
          
          
        }
      }
      
      Sigma_theta<-t(UU)%*%UU
      
      # Sigma_theta<-cov2cor(Sigma_theta)
    }
    
    # Interc_H<- Interc_HT
    #Slope_H[,1:D]<-Slope_HT[,1:D]
    
    
    eps<-max(abs( Pre_Sigma_theta- Sigma_theta),abs( Pre_LAMBDA-LAMBDA),
             abs( Pre_Beta- Beta))
    
    if(eps<epsCheck) passing=passing+1
    
    # Interc_H[1]<-Interc_HT[1]
    i<-i+1
    IncreSize=IncreSize+5
    
  }
  
  Slope_B_Pre<-Slope_B
  
  # stopCluster(cl)
  
  
  
  
  
  
  if( ReturnLikelihood){
    AS<-t(Slope_H)%*%Sigma_theta%*%Slope_H
    
    SdAS<-sqrt(1+diag(AS))
    SdASM[,]<-I1L%*%t(SdAS)
    TsAS<-SdAS%*%t(SdAS)
    # TsASA[,,]<-TsAS
    
    AT<-AS/TsAS
    diag( AT)<-1 #cor of diagonal equal to 1
    
    ST[,,]<- AT
    #SA is covariance matrix
    SA<-IA2*ST
    
    #TsASA[,,]<-TsAS
    # ASA[,,]<-AS
    # SA<-(ASA*IA2)/TsASA
    
    #for(l in 1:L){
    #  diag(SA[,,l])<-1
    #}
    
    IAd<-IA*(I1L%*%t(Interc_H))
    Fq<-IAd/SdASM
    
    # for(l in 1:L){
    #   PA[l]<- 
    #     mvtnorm::pmvnorm(lower=LOW,upper=Fq[l,],mean=MEAN,corr=SA[,,l])
    # }
    
    
    
    
    PA<-apply(matrix(1:L),1,FUN=PACDF,LOW=LOW,upperMat=Fq,corrArr=SA)  
    
    
    
    PAM[,]<-PA%*%t(I1N)
    
    
    Interc_B_matrix[,]<-I1L%*%t(Interc_B)
    t_Beta1_Mat<-t(Slope_B)
    
    AB<-Poss_AttrUp%*%t_Beta1_Mat+Interc_B_matrix
    #########################
    #  AB[AB<0]<- 0
    ##########################
    
    PR1<-plogis(Poss_AttrUp%*%t_Beta1_Mat+Interc_B_matrix)
    
    PRt1<- PR1[rLN,]
    PRt2<-1-PRt1
    PRT1<-ResM*PRt1+ResMN*PRt2
    
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
    
    
    
    NUM<-PRTNJ*PAM
    DEM<-colSums(NUM)
    LIKELI<- sum(log(DEM))
    
    
  }else{
    LIKELI<-NULL
  }
  
  time_end<-proc.time()-time_start
  return(list(Slope_H=Slope_H,Interc_H=Interc_H,Sigma_theta=Sigma_theta,
              Interc_B=Interc_B,Slope_B=Slope_B,
              Slope_B_Pre=Slope_B_Pre,PRAA=PRAA,
              time_end=time_end,i=i,LIKELI=LIKELI))
  
  
}
