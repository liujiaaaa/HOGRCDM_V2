LoadingPackages<-function(){
  PACKAGES<-c("parallel", "stats","graphics","grDevices","utils","datasets",
              "methods","base","cluster","lpSolve","DescTools",
              "mvtnorm", "RcppHungarian","glmnet","Matrix",
              "doParallel", "iterators","foreach",
              "profvis","emdbook","mvnfast","matrixStats", "clue",
              "MASS","remotes","future.apply","magrittr", "TruncatedNormal"
  )
  
  
  for(lu in 1:length(PACKAGES)){
    if( !(PACKAGES[lu] %in% installed.packages())) {
      install.packages(PACKAGES[lu],lib = Sys.getenv("R_LIBS_USER"),
                       repos="https://mirror.las.iastate.edu/CRAN/")
      cat(paste0("Installing Package:", PACKAGES[lu]))
    }
  }
  
  for(lu in 1:length(PACKAGES)){
    library(PACKAGES[lu],
            character.only = TRUE, logical.return = TRUE)
  }
  
}
