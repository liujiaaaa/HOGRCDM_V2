# This function checks the input parameter for formulating model
ModelSetCheck <- function(ModelSetList) {
  
  ModelDistri_Bottom=ModelSetList$ModelDistri_Bottom
  ModelStruc_Bottom=ModelSetList$ModelStruc_Bottom
  ModelStruc_Higher=ModelSetList$ModelStruc_Higher
  ModelFrame_Bottom=ModelSetList$ModelFrame_Bottom
  ModelFrame_Higher=ModelSetList$ModelFrame_Higher
  HigherLayer=ModelSetList$ HigherLayer
  
  
  ModelDistri_BottomList <- c("Lognormal", "LLM", "Poisson","Gamma")
  ModelStruc_BottomList <- c("Main_effect", "All_effect", "DINA")
  ModelStruc_HigherList <- c("Subscale", "Bifactor")
  ModelFrame_BottomList<-c("Confirmatory", "Exploratory")
  ModelFrame_HigherList<-c("Confirmatory", "Exploratory")
  
  if(ModelSetList$HigherLayer==T){
    cat("Higher order CDM")
  }else{
    cat("Single layer CDM")
  }
  
  
  
  # Check if ModelDistri_Bottom is in ModelDistri_BottomList
  if (ModelDistri_Bottom %in% ModelDistri_BottomList) {
    print(paste("Bottom Layer Model Distribution:", ModelDistri_Bottom))
  } else {
    stop("Error: ModelDistri_Bottom must be one of ", paste(ModelDistri_BottomList, collapse = ", "))
  }
  
  
  if (ModelStruc_Bottom %in% ModelStruc_BottomList) {
    print(paste("Bottom Layer Model Structure:", ModelStruc_Bottom))
  } else {
    stop("Error: ModelStruc_Bottom must be one of ", paste(ModelStruc_BottomList, collapse = ", "))
  }
  
  
  if (ModelFrame_Bottom %in% ModelFrame_BottomList) {
    print(paste("Bottom Layer Model Framework:", ModelFrame_Bottom))
  } else {
    stop("Error: ModelFrame_Bottom must be one of ", paste( ModelFrame_BottomList, collapse = ", "))
  }
  
  if(ModelSetList$HigherLayer==T){
     if (ModelStruc_Higher %in% ModelStruc_HigherList) {
    print(paste("Higher Layer Model Structure:", ModelStruc_Higher))
  } else {
    stop("Error: ModelStruc_Higher must be one of ", paste(ModelStruc_HigherList, collapse = ", "))
  }
  
  if (ModelFrame_Higher %in% ModelFrame_HigherList) {
    print(paste("Higher Layer Model Framework:", ModelFrame_Higher))
  } else {
    stop("Error: ModelFrame_Higher must be one of ", paste( ModelFrame_HigherList, collapse = ", "))
  }
  }
  
 
  
  
  

  return(list(ModelDistri_Bottom=ModelDistri_Bottom,ModelStruc_Bottom=ModelStruc_Bottom,ModelStruc_Higher=ModelStruc_Higher,
              ModelFrame_Bottom= ModelFrame_Bottom, ModelFrame_Higher= ModelFrame_Higher,HigherLayer=HigherLayer))
  
}
