OptBott_Confir_speed <- function(j, x, Res, Q_mat, weights, Family, Link) {
  
  # Select the relevant columns based on Q_mat
  data_x <- x[, which(Q_mat[j,] != 0)]
  data_x_with_intercept <- cbind(1, data_x)  # Add a column of ones for the intercept
 
  
  MDATA<- Merge_Data(data_y=Res,data_x=data_x_with_intercept, weights)
 UpRes<-MDATA$data_y
 data_x_with_intercept<-MDATA$data_x
 weights<-MDATA$weights
  
  # Dynamically retrieve the family function
  family_function <- match.fun(Family)
  
  # Dynamically set the link function if specified
  family_with_link <- family_function(link = Link)
  
  # Fit the model using speedglm.wfit
  fit <- speedglm.wfit(y =  UpRes[, j], X = data_x_with_intercept, family = family_with_link,
                       weights = weights, intercept = TRUE)
  
  # Extract coefficients
  res <- as.vector(coef(fit))
  
  # Prepare the result vector
  VecF <- rep(0, ncol(Q_mat) + 1)
  VecF[1] <- res[1]  # Intercept
  VecF[-1][which(Q_mat[j,] != 0)] <- res[-1]  # Coefficients
  
  # Return the result
  VecF
}

