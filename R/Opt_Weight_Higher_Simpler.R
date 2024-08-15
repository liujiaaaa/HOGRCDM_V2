Opt_Weight_Higher_Simpler<-function(j,x,Q_mat,Res,Link="probit", weights){
  data_x <- x[, which(Q_mat[j,] != 0)]
data_x_with_intercept <- cbind(1, data_x) # Add a column of ones for the intercept
 #data_x_with_intercept<-data_x
  
  fit <- glm.fit(x = data_x_with_intercept, y = Res[, j], weights = weights,intercept=T,
                 family = binomial(link = Link))
  
  res<-as.vector(coef(fit))
  
  VecF<-rep(0, ncol(Q_mat)+1)
  VecF[1]<-res[1]
  VecF[-1][which(Q_mat[j,] != 0)]<-res[-1]
  VecF
}

# test<-  glm.fit(x=thetaD, y=xxD[,1], weights = wwe,
#                 family = binomial(link = "probit"),
#                 control = list(), intercept = TRUE, singular.ok = TRUE)