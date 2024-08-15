OptBott_LogNor_Confir<-function(j,x,Res,Q_mat, weights){
  
  data_x <- x[, which(Q_mat[j,] != 0)]
  data_x_with_intercept <- cbind(1, data_x) # Add a column of ones for the intercept
  
  
  fit <- glm.fit(x = data_x_with_intercept, y = Res[, j], weights = weights,intercept=T,
                 family = gaussian())
  res<-as.vector(coef(fit))
  
 VecF<-rep(0, ncol(Q_mat)+1)
 VecF[1]<-res[1]
 VecF[-1][which(Q_mat[j,] != 0)]<-res[-1]
 VecF
  #parallel=T
  
  
  
}
# j=1
# fit <- glmnet(x=Data_Attri0, y=Data_Layer0[,j],  weights=WEI1,family="binomial",
#              )
# as.vector(coef(fit))

