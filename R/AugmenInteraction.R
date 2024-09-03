
#' Augmentation
#'
#' Expands the original binary latent attribute matrix to include all interaction 
#' terms, which is required for the "All_effect" and "DINA" cases.
#' 
#' @param OrigriData An original binary attribute matrix with K columns.
#' @return A list containing the augmented matrix 'OrigriData_Aug' and the number 
#' of columns in 'OrigriData_Aug', 's_num'.
#' @examples 
#' A_test <- matrix(c(1,1,1,0,0,1,0,0), ncol=2, byrow=TRUE)
#' A_test_Aug <- AugmenInteraction(A_test, 2)
#' @export






AugmenInteraction<-function(OrigriData,K){

  OrigriData_Aug<-OrigriData
  s_num<-K
  for(kk in 2:K){
    combinations <- t(combn(1:K, kk))
    te<-apply(OrigriData,1,FUN= comFUN, combinations)
    # print(dim(as.matrix(te)))
    # print(choose(K,kk))
    if(kk<K) te<-t(te)
    OrigriData_Aug<-cbind(OrigriData_Aug,te)
    rownames( OrigriData_Aug)<-NULL
    colnames( OrigriData_Aug)<-NULL
    s_num<-s_num+choose(K,kk)
  }
  return(list(OrigriData_Aug=OrigriData_Aug,s_num=s_num))
}

