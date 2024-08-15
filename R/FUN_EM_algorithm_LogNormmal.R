
FUN_EM_algorithm_LogNormmal<-function(Res, size_list, Initial_List, SettingList,gamma, Q_Higher,
                                      burnin=10,LLambda,Com_par=NULL,isBifactor=F,
                                      SubSampleEStep=F,Pre_period=10,Subsize=1,
                                      Likelihood=F,EXPECTATION=F,PACKAGES=NULL,type="ACDM",ConP=NULL){
  
  
 
  
  # install.packages("TruncatedNormal",lib = custom_library_path,
  #                  repos="https://mirror.las.iastate.edu/CRAN/")
 
  
  # # Load the package from the custom library pathz
  # library(remotes, lib.loc = custom_library_path)
  # #library(remotes, lib.loc = custom_library_path)
  # 
  # 
  # #  remotes::install_github("DavisVaughan/furrr")
  # library(furrr,lib.loc = custom_library_path)
  
  
  Opt_Weight_Higher<-function(j,x,Q_mat,Res,LLambda,Link="identity", weights){
    data_x<-x[,which(Q_mat[j,]!=0)]
    fit <- glmnet(x=data_x,Res[,j],  family = binomial(link = Link), weights=weights,
                  lambda=LLambda[j])
    as.vector(coef(fit))
  }
  
  
  
  
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
  
  
  
  
  
  
  
  sample_thetas <- function(L,K,D,data,Lam_1,Lam_0,Sigma_theta,  mc_samples){
    #source("sample_single_thetas.R")
    # source("/burg/home/jl6795/DINA2/DINA/sample_single_thetas.R")
    n <- L
    m <- K
    k <- D
    library(magrittr)
    
    d1_array <- sweep((rep(1, n) %x% t( Lam_1)), 1, c(t(2*data-1)), "*")
    d2_array <- sweep(2*data-1, 2, (Lam_0), "*")
    
    s_array<- sqrt(diag(t(Lam_1)%*%Sigma_theta%*%Lam_1)+1)
    
    # s_array <- (rowSums(t(Lam_1)^2)+1)^0.5
    thetas <- 1:n %>% 
      furrr::future_map(function(.x) sample_single_thetas(m, k, mc_samples, d1_array[((.x-1)*m+1):(.x*m), ], d2_array[.x, ], s_array,Sigma_theta), .options = furrr::furrr_options(seed = TRUE)) %>% 
      do.call(rbind, .)
    return(thetas)
  }
  
  
  time_strat<-proc.time()
  # library(mvnfast)
  # library(magrittr)
  # install.packages("remotes",lib = Sys.getenv("R_LIBS_USER"),
  #                  repos="https://mirror.las.iastate.edu/CRAN/")
  # install.packages("furrr",lib = Sys.getenv("R_LIBS_USER"),
  #                  repos="https://mirror.las.iastate.edu/CRAN/")
  
  
  
  N<-size_list$N
  J<-size_list$J
  K<-size_list$K
  D<-size_list$D
  maxIter<-size_list$maxIter
  L<-2^K
  lRes<-log(Res)
  
  
  Alpha<-Initial_List$Alpha0
  #theta<-Initial_List$theta0
  Lam_1<-Initial_List$Lam_1_0
  Lam_0<-Initial_List$Lam_0_0
  Beta0<-Initial_List$Beta0_0
  Beta1_Mat<-Initial_List$Beta1_Mat0
  
  if(is.null(ConP)){
    ConP=SettingList$ConP
  }
  
  
  if(type=="GDINA"|type=="DINA"){
    Q_Aug<-size_list$Q_Aug
    InitListUpdate<-list()
    InitListUpdate$Beta1_Mat<-matrix(0,nrow(Q_Aug),ncol(Q_Aug))
    InitListUpdate$Beta1_Mat[,1:K]<-Beta1_Mat
    InitListUpdate$Beta0<-Beta0
    
    
    Beta0<-InitListUpdate$Beta0
    Beta1_Mat<-InitListUpdate$Beta1_Mat
  }
  
  Sigma_theta<-Initial_List$Sigma_theta0
  ta<-Initial_List$sdd0#initial values for sd
  
  # Lam_1<- Lam_1T
  # Lam_0<- Lam_0T
  
  #Beta1_Mat<-Beta1_Mat*Q
  
  isBifactor<-SettingList$isBifactor
  Pre_period<-SettingList$Pre_period
  Subsize<-SettingList$Subsize
  maxIter<-SettingList$maxIter
  burnin<-SettingList$burnin
  gamma<-SettingList$gamma
  LLambda<-SettingList$LLambda
  Pre_period2<-SettingList$Pre_period2
  Likelihood<-SettingList$Likelihood
  
  
  #Lam_1<-Lam_1T
  #Lam_0<-Lam_0T
  #  Beta1_Mat<-Beta1_MatT
  # Beta0<-Beta0T
  
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
  
  #Data_x<-Poss_Attr
  #for(su in 2:N){
  #  Data_x<-rbind(Data_x,Poss_Attr)
  #}
  if(type=="ACDM"){
    Data_x1<- matrix(NA,K,L*N)
    Data_x1[,]<-t(Poss_Attr)
    Data_x<-t(Data_x1)
    
    Data_Res<-matrix(NA,N*L,J)
    for(l in 1:L){
      Data_Res[ ((0:(N-1))*L)+l,]<-lRes
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
    
    Data_Res<-DataSet_Create(num_ob=N,num_attri=K,num_item=J,OriData=lRes)
    
    Poss_AttrUp<- Poss_Attr1_Aug
  }
  
  
  
  #Data_Res<-matrix(NA,N*2^K,J)
  #for(su in 1:N){
  #  Data_Res[((su-1)*L+1):(su*L),]<-rep(1,L)%*%t(as.vector(Res[su,]))
  #}
  

  
  
  Lambda<-rbind(Lam_1,Lam_0)
  Lambda[-(D+1),][t(Q_Higher)==0]<-0
  QQc<-cbind(Q_Higher,rep(1,K))
  
  #This is used for the parameter estimation 
  index<-list()
  for(k in 1:K){
    index[[k]]<- which(Q_Higher[k,]==1)
  }
  
  
  # Create the T4 list with zeros
  #  T4 <- replicate(J, 0, simplify = FALSE)
  # T5<-T6 <- replicate(K, 0, simplify = FALSE)
  
  T7<-0
  T1<-T2<-0
  
  #T4<-list()
  #for(j in 1:J){
  # T4[[j]]<-0
  #}
  
  
  
  
  
  #T5<-T6<-list()
  #for(k in 1:K){
  # T5[[k]]<-T6[[k]]<-0
  #}
  
  
  
  
  #IA2<-array(NA,dim=c(K,K,L))
  #for(l in 1:L){
  #  for(ii in 1:K){
  #    for(jj in 1:K){
  #      IA2[ii,jj,l]<-
  #       (-1)^(NPoss_Attr[l,ii]+NPoss_Attr[l,jj])
  #   }
  # }
  #}
  
  IA2<- array((-1)^t((NPoss_Attr[, rep(1:K, each = K)] +
                        NPoss_Attr[, rep(1:K, times = K)])), dim = c(K, K, L))
  
  
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
  # for(l in 1:L){
  #   ResA[l,,]<-t(Res)
  # }
  # ResAN<-1-ResA
  
  lResM<-lRes[rep(1:nrow(Res),L),]
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
  
  
  ###############Fix the initial burnin=5 and increase it by 5
  ###############The previously specified burnin won't work
  burnin=5
  ################
  
  Beta<-rbind( Beta0,t(Beta1_Mat) )
  passing=0
  
  while((passing<3&i<maxIter)){
    
    
    # if(i%%10==0){
    #      print(Beta0)
    # print( Beta1_Mat)
    # print( Lam_1)
    # print( Lam_0)
    # print(Sigma_theta)
    # print(ta)
    # }
    


    Pre_Lambda<-Lambda
    Pre_Sigma_theta<-Sigma_theta
    Pre_Beta<-  Beta
    
    AS<-t(Lam_1)%*%Sigma_theta%*%Lam_1
    
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
    
    IAd<-IA*(I1L%*%t(Lam_0))
    Fq<-IAd/SdASM
    
    # for(l in 1:L){
    #   PA[l]<- mvtnorm::pmvnorm(lower=LOW,upper=Fq[l,],mean=MEAN,corr=SA[,,l])
    # }
    
   
    
  PA<-apply(matrix(1:L),1,FUN=PACDF,LOW=LOW,upperMat=Fq,corrArr=SA)  
    
  ###################################################  
   if(i<ConP)  PA<-PA*1/2+rep(1/L,L)*1/2
  ###################################################  
    
    
    
    #  PA <- foreach(l = 1:L, .combine = "c") %dopar% {
    #   library(mvtnorm)
    #   pmvnorm(lower = LOW, upper = Fq[l,], mean = MEAN, corr = SA[,,l])
    # }
    
    
    PAM[,]<-PA%*%t(I1N)
    
    
    Beta0_matrix[,]<-I1L%*%t(Beta0)
    t_Beta1_Mat<-t(Beta1_Mat)
    
    AB<-Poss_AttrUp%*%t_Beta1_Mat+Beta0_matrix
    #########################
    #  AB[AB<0]<- 0
    ##########################
    
    sdd<-rep(ta,each=N*L)
    PRT1[,]<-dnorm(as.vector(lResM),as.vector(AB[rLN,]),sd= sdd)
    
    
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
    DEMM<-rep(1,L)%*%t(DEM)
    
    #Posterior
    PRAA<-NUM/DEMM 
    #default family is gaussian
    
    
    
    
    
    ppA<-(rowSums(PRAA))
    
    PRAV<-sqrt(ppA)
    PRAVM<-PRAV%*%t(rep(1,D))
    WEI<-as.vector(PRAA)
    
    PRAV_s<-as.vector(t(PRAA))%*%t(rep(1,J))
    ta<- sqrt( colSums((lResM-AB[rLN,])^2* PRAV_s)/N)
    
    
    
    
   
    
    
    
    Beta<-sapply(1:J,Opt_WeightFast, x=Data_x, 
                 Res=Data_Res,weights=WEI, 
                 LLambda=LLambda, simplify = T)
    
    
    
    Beta0 <- Beta[1,]
    Beta1_Mat <- t(Beta[-1,])
    
    
    
    colnames(Beta0)<-NULL
    rownames(Beta1_Mat)<-colnames(Beta1_Mat)<-NULL
    
    
  
    
    thetas1<- sample_thetas(L,K,D,data=Poss_Attr,Lam_1,Lam_0,
                            Sigma_theta=Sigma_theta,mc_samples=burnin)
    thetaA<-array(NA,dim=c(L,burnin,D))
    for(ll in 1:L){
      # print((((ll-1)*mc_samples)+1):(ll*mc_samples))
      thetaA[ll,,]<-thetas1[(((ll-1)*burnin)+1):(ll*burnin),]
    }
    thetaA1<-aperm( thetaA,c(3,2,1))
    
  
    theta_At<-matrix(NA,D,L*burnin)
    theta_At[,]<-thetaA1
    thetaD<-t(theta_At[,])
    
    xxD<-Poss_Attr[rep(1:L,each=burnin),]
    
    wwe<-rep(ppA,each=burnin)   
    
    wee1<- rep(PRAV,each=burnin)  
    thetas2<- thetas1*(wee1%*%t(rep(1,D)))
  
    
    
    
    
    # LAMBDA<-sapply(1:K,  Opt_Weight_Higher,Q_mat=Q_Higher,
    #                x=thetaD,
    #                Res=xxD,weights= wwe, Link="probit",
    #                LLambda=rep(0,K), simplify = T)
    
    LAMBDA<-sapply(1:K,    Opt_Weight_Higher_Simpler,Q_mat=Q_Higher,
                   x=thetaD,
                   Res=xxD,weights= wwe, Link="probit",
                   LLambda=rep(0,K), simplify = T)
    

    
    #print(test)
    
    result_Lamnda<-rbind(LAMBDA[-1,],LAMBDA[1,])
    Lambda[t(QQc)!=0]<-result_Lamnda
  #  print( Lambda)
    
    
    Lam_1<-Lambda[-(D+1),]
    Lam_0<-Lambda[(D+1),]  
    
   # Sigma_theta<-(theta_At)%*%t(theta_At)/(L*burnin)
    Sigma_theta<-t(thetas2)%*%thetas2/(burnin*N)
    
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
    
    # Lam_0<- Lam_0T
    #Lam_1[,1:D]<-Lam_1T[,1:D]
    
    
    eps<-max(abs( Pre_Sigma_theta- Sigma_theta),abs( Pre_Lambda-Lambda),
             abs( Pre_Beta- Beta))
    
    if(eps<0.04) passing=passing+1
    
    # Lam_0[1]<-Lam_0T[1]
    i<-i+1
    burnin=burnin+5
    
  }
  
  Beta1_Mat_Pre<-Beta1_Mat
  
  #stopCluster(cl)
  
  
  KK=ncol(Com_par$Beta1_MatT)
  DD<-nrow(Com_par$Lam_1T)
  if(!is.null(Com_par)){
    
    Cost_matrix<-matrix(NA,KK,KK)
    for(ll in 1:KK){
      for(mm in 1:KK){
        Cost_matrix[ll,mm]<-sum(abs(Beta1_Mat[,mm]-Com_par$Beta1_MatT[,ll]))
      } 
    }
    Test<-HungarianSolver(Cost_matrix) 
    if(sum(abs(Test$pairs[,2]-(1:KK)))!=0) {
      print(Test$pairs[,2]);Hungarian=T
    }else{
      Hungarian=F
    }
    Beta1_Mat<-Beta1_Mat[, Test$pairs[,2]]
    # Lam_1<-Lam_1[,Test$pairs[,2]]
    # Lam_0<-Lam_0[Test$pairs[,2]]
    # 
    # 
    # Cost_matrixH<-matrix(NA,DD,DD)
    # for(ll in 1:DD){
    #   for(mm in 1:DD){
    #     Cost_matrixH[ll,mm]<-sum(abs(Lam_1[mm,]-Com_par$Lam_1T[ll,]))
    #   }
    # }
    # 
    # TestH<-HungarianSolver(Cost_matrixH)
    # if(sum(abs(TestH$pairs[,2]-(1:DD)))!=0) {
    #   print(TestH$pairs[,2]);Hungarian=T
    # }else{
    #   Hungarian=F
    # }
    # Lam_1<-Lam_1[TestH$pairs[,2],]
    
    
    
  }
 
  
  
  if(Likelihood){
    AS<-t(Lam_1)%*%Sigma_theta%*%Lam_1
    
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
    
    IAd<-IA*(I1L%*%t(Lam_0))
    Fq<-IAd/SdASM
    
    # for(l in 1:L){
    #   PA[l]<- 
    #     mvtnorm::pmvnorm(lower=LOW,upper=Fq[l,],mean=MEAN,corr=SA[,,l])
    # }
    

    
    
    PA<-apply(matrix(1:L),1,FUN=PACDF,LOW=LOW,upperMat=Fq,corrArr=SA)  
    
    
    
    PAM[,]<-PA%*%t(I1N)
    
    
    Beta0_matrix[,]<-I1L%*%t(Beta0)
    t_Beta1_Mat<-t(Beta1_Mat)
    
    AB<-Poss_AttrUp%*%t_Beta1_Mat+Beta0_matrix
    #########################
    # AB[AB<0]<- 0
    ##########################
    
    PRT1[,]<-dnorm(as.vector(lResM),as.vector(AB[rLN,]),sdd)
    
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
  
  time_end<-proc.time()-time_strat
  return(list(Lam_1=Lam_1,Lam_0=Lam_0,Sigma_theta=Sigma_theta,
              Beta0=Beta0,Beta1_Mat=Beta1_Mat,
              Beta1_Mat_Pre=Beta1_Mat_Pre, Pairs=Test$pairs,PRAA=PRAA,
              Hungarian=Hungarian,time_end=time_end,i=i,LIKELI=LIKELI,SD= ta))
  
  
}
