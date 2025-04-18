OptBott_Confir <- function(j, x, Res, Q_mat, weights, Family, Link) {
  
  # Select the relevant columns based on Q_mat
  data_x <- x[, which(Q_mat[j,] != 0)]
  data_x_with_intercept <- cbind(1, data_x) # Add a column of ones for the intercept
  
  # Dynamically retrieve the family function
  family_function <- match.fun(Family)
  
  # Dynamically set the link function if specified
  family_with_link <- family_function(link = Link)
  
  # Fit the model using glm.fit
  fit <- glm.fit(x = data_x_with_intercept, y = Res[, j], weights = weights, intercept = TRUE,
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
