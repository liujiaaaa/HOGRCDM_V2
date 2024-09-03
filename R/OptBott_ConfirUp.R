OptBott_ConfirUp <- function(j, DATA_INDEX_LIST,Q_mat, weights, N,Family, Link) {
  
  x=as.matrix(DATA_INDEX_LIST[[j]]$dX)
  Res1= as.matrix(DATA_INDEX_LIST[[j]]$dRes)
  weightUp<-tapply(weights,DATA_INDEX_LIST[[j]]$INDEX_Wei, sum)/N
  weightUp<-as.vector(weightUp)
  
  inR<-rep(1:nrow(Res1),ceiling(100/nrow(Res1)))
  inRUP<-sample(inR)
  
  
 # SI<-sample(1:nrow(Res1),size=20000, replace=T)
  
  # Select the relevant columns based on Q_mat
  data_x <-x[inRUP,]
  data_x_with_intercept <- cbind(1, data_x) # Add a column of ones for the intercept
 # weightUp1<-weightUp
   weightUp1<- as.vector(weightUp[ inRUP])
   weightUp1<-weightUp1/sum(weightUp1)                  
  
   ResUp<-as.matrix(Res1[inRUP])
  # Dynamically retrieve the family function
  family_function <- match.fun(Family)
  
  # Dynamically set the link function if specified
  family_with_link <- family_function(link = Link)
  
  # Fit the model using glm.fit
  fit <- glm.fit(x = data_x_with_intercept, y = ResUp,  intercept = TRUE, weights=weightUp1,
                 family = family_with_link)
  
  # Extract coefficients
  res <- as.vector(coef(fit))
  
  # Prepare the result vector
  VecF <- rep(0, ncol(Q_mat) + 1)
  VecF[1] <- res[1]
  VecF[-1][which(Q_mat[j,] != 0)] <- res[-1]
  
  # Return the result
  VecF
}

# Example usage:
# OptBott_LogNor_Confir(j, x, Res, Q_mat, weights, "gaussian", "identity")
# OptBott_LogNor_Confir(j, x, Res, Q_mat, weights, "binomial", "logit")


# Example usage:
# OptBott_LogNor_Confir(j, x, Res, Q_mat, weights, "gaussian")

# j=1
# fit <- glmnet(x=Data_Attri0, y=Data_Layer0[,j],  weights=WEI1,family="binomial",
#              )
# as.vector(coef(fit))
