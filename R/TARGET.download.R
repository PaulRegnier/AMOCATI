#' Downloads the data (RNA-Seq + clinical) from the desired TARGET project
#'
#' This function downloads the data (RNA-Seq + clinical) from the desired GDC TARGET project.
#'
#' @param projectID A character value determining which project to use.
#' 
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @return Generated output files are fist exported to the `temp` folder before being moved to the `output > data` folder.
#'
#' @importFrom foreach %do%
#' @importFrom dplyr %>%
#'
#' @export

TARGET.download = function(projectID = NULL, verbose = TRUE)
{
  if(verbose == TRUE)
  {
    print("## RNASeq data metamapping ##")
  }
  
  metaMappingRNASeqRaw = GenomicDataCommons::files() %>% GenomicDataCommons::filter(~ cases.project.project_id == projectID & data_type == "Gene Expression Quantification" & access == "open") %>% GenomicDataCommons::select(c("file_id", "cases.case_id", "data_type", "file_name", "file_size")) %>% GenomicDataCommons::response_all()
  
  metaMappingRNASeqManifest = metaMappingRNASeqRaw %>% GenomicDataCommons::manifest()
  
  metaMappingRNASeqRawRemoveIDIncrement = gsub("_id(.+)", "_id", names(unlist(metaMappingRNASeqRaw$results$cases)))
  RNASeq_IDsToKeep = which(duplicated(metaMappingRNASeqRawRemoveIDIncrement) == FALSE)
  metaMappingRNASeqRawNoDuplicates = as.vector(unlist(metaMappingRNASeqRaw$results$cases)[RNASeq_IDsToKeep])
  metaMappingRNASeqRaw$results$cases = metaMappingRNASeqRawNoDuplicates
  
  
  utils::write.table(metaMappingRNASeqRaw$results, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "metaMappingRNASeqRaw.txt"), sep = "\t")
  utils::write.table(metaMappingRNASeqManifest, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "metaMappingRNASeqRaw-manifest.txt"), sep = "\t")
  
  if(verbose == TRUE)
  {
    print("## Clinical data metamapping ##")
  }
  
  metaMappingClinicalRaw = GenomicDataCommons::files() %>% GenomicDataCommons::filter(~ cases.project.project_id == projectID & data_type == "Clinical Supplement" & access == "open") %>% GenomicDataCommons::select(c("file_id", "data_type", "file_name", "file_size")) %>% GenomicDataCommons::response_all()
  
  metaMappingClinicalRawManifest = metaMappingClinicalRaw %>% GenomicDataCommons::manifest()
  
  utils::write.table(metaMappingClinicalRaw$results, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "metaMappingClinicalRaw.txt"), sep = "\t")
  utils::write.table(metaMappingClinicalRawManifest, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "metaMappingClinicalRaw-manifest.txt"), sep = "\t")
  
  if(verbose == TRUE)
  {
    print("## Emptying final data folder ##")
  }
  
  unlink(file.path("data", "RNASeq", "*.*"), TRUE)
  unlink(file.path("data", "Clinical", "*.*"), TRUE)
  
  if(verbose == TRUE)
  {
    print("## RNASeq and clinical data downloading to cache folder ##")
  }
  
  tempFolder = file.path("temp")
  GenomicDataCommons::gdc_set_cache(tempFolder)
  
  if(verbose == TRUE)
  {
    print("## RNASeq data downloading to cache folder (may be long) ##")
  }
  
  downloadedRNASeqFiles = GenomicDataCommons::gdcdata(metaMappingRNASeqManifest$id, use_cached = TRUE, progress = TRUE)	
  downloadedRNASeqFilesList = names(downloadedRNASeqFiles)
  
  if(verbose == TRUE)
  {
    print("## Unzipping and moving RNASeq data files to destination folder ##")
  }
  
  a = NULL
  
  foreach::foreach(a = 1:length(downloadedRNASeqFilesList)) %do%
    {
      currentFolderToOpen = file.path(tempFolder, downloadedRNASeqFilesList[a])
      currentFileToMove = list.files(currentFolderToOpen)[1]
      currentPathToMove = file.path(currentFolderToOpen, currentFileToMove)
      currentPathDestination = file.path("data", "RNASeq", gsub(".gz", "", currentFileToMove))
      
      R.utils::gunzip(filename = currentPathToMove, destname = currentPathDestination, overwrite = TRUE, remove = FALSE)
      #unlink(currentFolderToOpen)
    }
  
  if(verbose == TRUE)
  {
    print("## Clinical data downloading to cache folder (may be long) ##")
  }
  
  downloadedClinicalFiles = GenomicDataCommons::gdcdata(metaMappingClinicalRawManifest$id, use_cached = TRUE, progress = TRUE)
  downloadedClinicalFilesList = names(downloadedClinicalFiles)
  
  if(verbose == TRUE)
  {
    print("## Moving clinical data files to destination folder ##")
  }
  
  foreach::foreach(a = 1:length(downloadedClinicalFilesList)) %do%
    {
      currentFolderToOpen = file.path(tempFolder, downloadedClinicalFilesList[a])
      currentFileToMove = list.files(currentFolderToOpen)[1]
      currentPathToMove = file.path(currentFolderToOpen, currentFileToMove)
      currentPathDestination = file.path("data", "Clinical", paste("TARGETClinicalData_", a, ".xlsx", sep = ""))
      
      file.copy(from = currentPathToMove, to = currentPathDestination, overwrite = TRUE)
      #unlink(currentFolderToOpen)
    }
  
  if(verbose == TRUE)
  {
    print("## WARNING : open each clinical data file (xslx format) and check if you need to pool them or not. Then save the file as 'TARGETClinicalData.xslx' ##")
  }
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
}