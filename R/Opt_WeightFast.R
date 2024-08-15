Opt_WeightFast<-function(j,x,Res,LLambda, weights, family="gaussian",lower.limits=-Inf,
                         stan = TRUE,thresh = 1e-07){
  
  fit <- glmnet(x, y=Res[,j],  weights=weights, lambda=LLambda[j],family=family,
                lower.limits=lower.limits,parallel=T)
  as.vector(coef(fit))
  
  #parallel=T
 
}
# j=1
# fit <- glmnet(x=Data_Attri0, y=Data_Layer0[,j],  weights=WEI1,family="binomial",
#              )
# as.vector(coef(fit))

