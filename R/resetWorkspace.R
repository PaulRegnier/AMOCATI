#' Reset the workspace
#'
#' This function resets the workspace that AMOCATI will use for its tasks. It can also erase the whole R memory.
#'
#' @param eraseEntireRMemory A boolean value indicating if the entire R memory should be erased (if set to `TRUE`) or not (if set to `FALSE`).
#'
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @export

resetWorkspace = function(eraseEntireRMemory = FALSE, verbose = TRUE)
{
  if(verbose == TRUE)
  {
    print("## Reset of the workspace: emptying and recreating folders ##")
  }
  
  foldersToErase = list.dirs(file.path("output"))
  unlink(foldersToErase, force = TRUE, recursive = TRUE)
  
  Sys.sleep(0.5)
  
  dir.create(file.path("output"))
  dir.create(file.path("output", "apply"))
  dir.create(file.path("output", "apply", "input"))
  dir.create(file.path("output", "apply", "customSignatures"))
  dir.create(file.path("output", "apply", "customSignatures", "fullTables"))
  dir.create(file.path("output", "apply", "customSignatures", "fullTables", "noDG"))
  dir.create(file.path("output", "apply", "customSignatures", "synthesis"))
  dir.create(file.path("output", "apply", "customSignatures", "CS"))
  dir.create(file.path("output", "apply", "customSignatures", "QS"))
  dir.create(file.path("output", "apply", "customSignatures", "CS VS QS"))
  dir.create(file.path("output", "data"))
  dir.create(file.path("output", "data", "input"))
  dir.create(file.path("output", "metaResults"))
  dir.create(file.path("output", "metaResults", "selectedGenes"))
  dir.create(file.path("output", "metaResults", "selectedGenes", "survivalTables"))
  dir.create(file.path("output", "metaResults", "selectedGenes", "survivalCurves"))
  dir.create(file.path("output", "signatures"))
  dir.create(file.path("output", "class"))
  dir.create(file.path("output", "class", "input"))
  
  foldersToErase = list.dirs(file.path("data"))
  unlink(foldersToErase, force = TRUE, recursive = TRUE)
  
  Sys.sleep(0.5)
  
  dir.create(file.path("data"))
  dir.create(file.path("data", "RNASeq"))
  dir.create(file.path("data", "Clinical"))	
  
  Sys.sleep(0.5)
  
  foldersToErase = list.dirs(file.path("temp"))
  unlink(foldersToErase, force = TRUE, recursive = TRUE)
  
  Sys.sleep(0.5)
  
  dir.create(file.path("temp"))
  
  if(verbose == TRUE)
  {
    print("## Purging computer RAM ##")
  }
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
  
  if(verbose == TRUE & eraseEntireRMemory == TRUE)
  {
    print("## Erasing entire R workspace and objects in memory ##")
  }
  
  if(eraseEntireRMemory == TRUE)
  {
    rm(list = ls(all.names = TRUE))
  }
}